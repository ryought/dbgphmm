//!
//! Neighbor copy numbers related functions
//!
//!
//!

use super::{CopyNums, MultiDbg};
use crate::graph::k_shortest::{k_shortest_cycle, k_shortest_simple_path};
use crate::graph::utils::split_node;

use itertools::Itertools;
use petgraph::graph::{DefaultIx, DiGraph, EdgeIndex, NodeIndex};
use petgraph_algos::common::is_edge_simple;
use rustflow::min_flow::{
    base::{FlowEdgeBase, FlowGraph},
    enumerate_neighboring_flows, find_neighboring_flow_by_edge_change,
    residue::{
        flow_to_residue_convex, is_meaningful_move_on_residue_graph, residue_graph_cycle_to_flow,
        update_cycle_from_str, update_cycle_to_string, ResidueDirection, UpdateCycle,
    },
};

///
/// Parameters for neighbor search of copy numbers
///
#[derive(Clone, Debug, Copy)]
pub struct NeighborConfig {
    /// Max size of cycle (in compact graph) in BFS short-cycle search
    ///
    pub max_cycle_size: usize,
    /// Max number of flips (+/- or -/+ changes) in cycles on compact residue graph
    ///
    pub max_flip: usize,
    /// Augment short cycles with long cycles causing 0x -> 1x change
    ///
    pub use_long_cycles: bool,
    /// Ignore cyclic paths passing through the terminal node NNN
    /// because this changes the number of haplotypes.
    ///
    pub ignore_cycles_passing_terminal: bool,
    ///
    /// Long cycle with Down direction only to reduce high copy number edges
    ///
    pub use_reducers: bool,
}

///
/// UpdateInfo
///
#[derive(Clone, Debug, PartialEq)]
pub struct UpdateInfo {
    pub cycles: Vec<UpdateCycle>,
    pub method: UpdateMethod,
}

impl UpdateInfo {
    pub fn new(cycles: Vec<UpdateCycle>, method: UpdateMethod) -> UpdateInfo {
        assert!(cycles.len() == 1 || method == UpdateMethod::MultiMove);
        UpdateInfo { cycles, method }
    }
    pub fn cycle(&self) -> &UpdateCycle {
        assert!(self.method != UpdateMethod::MultiMove);
        &self.cycles[0]
    }
}

///
///
///
#[derive(Clone, Debug, PartialEq, Copy)]
pub enum UpdateMethod {
    Rescue,
    MultiMove,
    Short,
    Long,
    Reducer,
    Manual,
}

///
/// Display/FromStr of UpdateInfo
///
///
/// ## Requirements
///
/// * No spaces (because INSPECT file uses TSV format)
/// * No commas (because comma separates elements of Vec<UpdateInfo>)
///
impl std::fmt::Display for UpdateInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.method {
            UpdateMethod::Rescue => {
                write!(f, "R({})", update_cycle_to_string(&self.cycles[0]))
            }
            UpdateMethod::MultiMove => {
                write!(
                    f,
                    "M({})",
                    self.cycles
                        .iter()
                        .map(|cycle| update_cycle_to_string(cycle))
                        .join(":")
                )
            }
            _ => {
                let method = match self.method {
                    UpdateMethod::Short => 'S',
                    UpdateMethod::Long => 'L',
                    UpdateMethod::Reducer => 'E',
                    UpdateMethod::Manual => 'X',
                    _ => unreachable!(),
                };
                write!(f, "{}({})", method, update_cycle_to_string(&self.cycles[0]))
            }
        }
    }
}
impl std::str::FromStr for UpdateInfo {
    type Err = ();
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let first = text.chars().nth(0).ok_or(())?;
        let method = match first {
            'R' => UpdateMethod::Rescue,
            'M' => UpdateMethod::MultiMove,
            'S' => UpdateMethod::Short,
            'L' => UpdateMethod::Long,
            'E' => UpdateMethod::Reducer,
            'X' => UpdateMethod::Manual,
            _ => return Err(()),
        };
        let (_, s) = text.split_at(1);
        let s = s.strip_prefix('(').ok_or(())?;
        let s = s.strip_suffix(')').ok_or(())?;
        if method == UpdateMethod::MultiMove {
            let cycles: Option<Vec<UpdateCycle>> =
                s.split(':').map(|t| update_cycle_from_str(t)).collect();
            let cycles = cycles.ok_or(())?;
            Ok(UpdateInfo::new(cycles, method))
        } else {
            let cycle = update_cycle_from_str(s).ok_or(())?;
            Ok(UpdateInfo::new(vec![cycle], method))
        }
    }
}

impl MultiDbg {
    /// Get neighboring copy numbers
    ///
    ///
    pub fn to_neighbor_copy_nums_and_infos(
        &self,
        config: NeighborConfig,
    ) -> Vec<(CopyNums, UpdateInfo)> {
        // (1) enumerate all short cycles
        let mut short_cycles = self.to_short_neighbors(config.max_cycle_size, config.max_flip);
        eprintln!("short_cycles: {}", short_cycles.len());

        // (2) add long cycles for 0x -> 1x
        if config.use_long_cycles {
            let mut long_cycles = self.to_long_neighbors();
            eprintln!("long_cycles: {}", long_cycles.len());
            short_cycles.append(&mut long_cycles);
        }

        // (3) add long cycles for reducing high-copy edges
        if config.use_reducers {
            let mut reducers = self.to_reducer_neighbors();
            eprintln!("reducers: {}", reducers.len());
            short_cycles.append(&mut reducers);
        }

        short_cycles
    }
    ///
    /// Create flow network `FlowGraph` for enumerating neighboring copy numbers (flow)
    ///
    pub fn to_flow_network(&self) -> FlowGraph<usize> {
        // create flow network for rustflow
        let network = self.graph_compact().map(
            |_, _| (),
            |edge, _| {
                let copy_num = self.copy_num_of_edge_in_compact(edge);
                FlowEdgeBase::new(copy_num.saturating_sub(1), copy_num.saturating_add(1), 0.0)
            },
        );
        network
    }
    /// Rescue neighbors
    /// 0x
    pub fn to_rescue_neighbors(
        &self,
        k_non_zero: usize,
        k_zero: usize,
        weighted_by_copy_num: bool,
    ) -> Vec<(CopyNums, UpdateInfo)> {
        let mut ret = vec![];

        let mut n_zero = 0;
        let mut n_non_zero = 0;

        for (edge, _, _, _) in self.edges_compact() {
            if self.copy_num_of_edge_in_compact(edge) == 0 {
                let cycles_non_zero =
                    self.to_rescue_neighbors_for_edge(edge, k_non_zero, true, weighted_by_copy_num);
                if cycles_non_zero.is_empty() {
                    let cycles_zero = self.to_rescue_neighbors_for_edge(
                        edge,
                        k_zero,
                        false,
                        weighted_by_copy_num,
                    );
                    // eprintln!("e{} zero {}", edge.index(), cycles_zero.len());
                    n_zero += 1;
                    ret.extend_from_slice(&cycles_zero);
                } else {
                    ret.extend_from_slice(&cycles_non_zero);
                    // eprintln!("e{} nonzero {}", edge.index(), cycles_non_zero.len());
                    n_non_zero += 1;
                }
            }
        }

        eprintln!("rescue n_zero={} n_non_zero={}", n_zero, n_non_zero);

        ret
    }
    /// Rescue neighbors
    ///
    /// 0x
    pub fn to_rescue_neighbors_for_edge(
        &self,
        edge: EdgeIndex,
        k: usize,
        not_make_new_zero_edge: bool,
        weighted_by_copy_num: bool,
    ) -> Vec<(CopyNums, UpdateInfo)> {
        let copy_nums = self.get_copy_nums();
        let network = self.graph_compact().map(
            |_, _| (),
            |edge, _| {
                let copy_num = self.copy_num_of_edge_in_compact(edge);
                if copy_num == 0 {
                    FlowEdgeBase::new(0, copy_num.saturating_add(1), 0.0)
                } else if not_make_new_zero_edge {
                    FlowEdgeBase::new(1, copy_num.saturating_add(1), 0.0)
                } else {
                    FlowEdgeBase::new(copy_num.saturating_sub(1), copy_num.saturating_add(1), 0.0)
                }
            },
        );
        // println!("{:?}", petgraph::dot::Dot::with_config(&network, &[]));
        let mut rg = flow_to_residue_convex(&network, &copy_nums);
        //
        // split terminal node into two nodes to prohibit cycles passing the terminal node.
        if let Some(terminal) = self.terminal_node_compact() {
            split_node(&mut rg, terminal, None);
        }
        let edge_in_rg = rg
            .edge_indices()
            .find(|&e| rg[e].target == edge && rg[e].direction == ResidueDirection::Up)
            .unwrap();
        let (v, w) = rg.edge_endpoints(edge_in_rg).unwrap();
        let cycles = k_shortest_simple_path(&rg, w, v, k, |e| {
            if e == edge_in_rg {
                usize::MAX
            } else {
                if weighted_by_copy_num {
                    self.n_bases(rg[e].target)
                        / self.copy_num_of_edge_in_compact(rg[e].target).max(1)
                } else {
                    self.n_bases(rg[e].target)
                }
            }
        });

        cycles
            .into_iter()
            .map(|mut cycle| {
                cycle.insert(0, edge_in_rg);
                cycle
            })
            .filter(|cycle| is_edge_simple(&rg, &cycle))
            .map(|cycle| residue_graph_cycle_to_flow(&copy_nums, &rg, &cycle))
            .map(|(copy_nums, cycle)| {
                (
                    copy_nums,
                    UpdateInfo::new(vec![cycle], UpdateMethod::Rescue),
                )
            })
            .collect()
    }
    ///
    ///
    ///
    pub fn to_short_neighbors(
        &self,
        max_cycle_size: usize,
        max_flip: usize,
    ) -> Vec<(CopyNums, UpdateInfo)> {
        let network = self.to_flow_network();
        let copy_nums = self.get_copy_nums();

        enumerate_neighboring_flows(&network, &copy_nums, Some(max_cycle_size), Some(max_flip))
            .into_iter()
            .map(|(copy_nums, cycle)| {
                (copy_nums, UpdateInfo::new(vec![cycle], UpdateMethod::Short))
            })
            .collect()
    }
    ///
    ///
    ///
    pub fn to_long_neighbors(&self) -> Vec<(CopyNums, UpdateInfo)> {
        let network = self.to_flow_network();
        let copy_nums = self.get_copy_nums();

        // let max_copy_num = self.max_copy_num();
        self.graph_compact()
            .edge_indices()
            .filter(|&e| self.copy_num_of_edge_in_compact(e) == 0)
            .filter_map(|e| {
                find_neighboring_flow_by_edge_change(
                    &network,
                    &copy_nums,
                    e,
                    ResidueDirection::Up,
                    // weight function
                    |e| self.n_bases(e) / (self.copy_num_of_edge_in_compact(e) + 1),
                    // |e| max_copy_num - self.copy_num_of_edge_in_compact(e),
                )
            })
            // .filter(|(_copy_nums, info)| info.len() > config.max_cycle_size)
            .filter(|(_copy_nums, cycle)| !self.is_passing_terminal(&cycle))
            .map(|(copy_nums, cycle)| (copy_nums, UpdateInfo::new(vec![cycle], UpdateMethod::Long)))
            .collect()
    }
    pub fn to_reducer_neighbors(&self) -> Vec<(CopyNums, UpdateInfo)> {
        let network = self.graph_compact().map(
            |_, _| (),
            |edge, _| {
                let copy_num = self.copy_num_of_edge_in_compact(edge);
                if copy_num > 2 {
                    // reduceable
                    FlowEdgeBase::new(copy_num.saturating_sub(1), copy_num, 0.0)
                } else {
                    FlowEdgeBase::new(copy_num, copy_num, 0.0)
                }
            },
        );
        let copy_nums = self.get_copy_nums();
        enumerate_neighboring_flows(&network, &copy_nums, Some(100), Some(0))
            .into_iter()
            .map(|(copy_nums, cycle)| {
                (
                    copy_nums,
                    UpdateInfo::new(vec![cycle], UpdateMethod::Reducer),
                )
            })
            .collect()
    }
    pub fn is_passing_terminal(&self, cycle: &UpdateCycle) -> bool {
        cycle
            .iter()
            .any(|&(edge, _)| self.is_start_or_end_edge_compact(edge))
    }
    pub fn has_zero_to_one_change(&self, cycle: &UpdateCycle) -> bool {
        cycle.iter().any(|&(edge, dir)| {
            self.copy_num_of_edge_in_compact(edge) == 0 && dir == ResidueDirection::Up
        })
    }
    ///
    /// check if all edges in `next_update` do not appears in `updates`
    /// to avoid conflicts of the simultaneous updates.
    ///
    pub fn is_independent_update(
        &self,
        copy_nums: &CopyNums,
        cycles: &[UpdateCycle],
        next_cycle: &UpdateCycle,
    ) -> bool {
        for cycle in cycles {
            for (edge, _) in cycle {
                for (next_edge, _) in next_cycle {
                    if edge == next_edge {
                        return false;
                    }
                }
            }
        }
        true
    }
    ///
    ///
    ///
    pub fn apply_update_info_to_copy_nums(&self, copy_nums: &mut CopyNums, cycle: &UpdateCycle) {
        for &(edge, direction) in cycle {
            match direction {
                ResidueDirection::Up => {
                    copy_nums[edge] += 1;
                }
                ResidueDirection::Down => {
                    copy_nums[edge] -= 1;
                }
            }
        }
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::super::toy;
    use super::*;
    use crate::common::ei;

    #[test]
    fn update_info_serde() {
        let ui1 = UpdateInfo::new(
            vec![vec![
                (ei(0), ResidueDirection::Up),
                (ei(1), ResidueDirection::Up),
                (ei(2), ResidueDirection::Down),
            ]],
            UpdateMethod::Rescue,
        );
        let ui2 = UpdateInfo::new(
            vec![
                vec![
                    (ei(0), ResidueDirection::Up),
                    (ei(1), ResidueDirection::Up),
                    (ei(2), ResidueDirection::Down),
                ],
                vec![
                    (ei(10), ResidueDirection::Down),
                    (ei(11), ResidueDirection::Down),
                ],
            ],
            UpdateMethod::MultiMove,
        );
        println!("{}", ui1);
        println!("{}", ui2);

        assert_eq!(ui1, ui1.to_string().parse().unwrap());
        assert_eq!(ui2, ui2.to_string().parse().unwrap());
    }

    #[test]
    fn neighbors_for_toy() {
        let mut dbg = toy::intersection();
        dbg.show_graph_with_kmer();
        let neighbors = dbg.to_neighbor_copy_nums_and_infos(NeighborConfig {
            max_cycle_size: 10,
            max_flip: 0,
            use_long_cycles: false,
            ignore_cycles_passing_terminal: false,
            use_reducers: false,
        });
        for (copy_nums, update_info) in neighbors {
            println!("{} {:?}", copy_nums, update_info);
        }

        {
            let neighbors = dbg.to_neighbor_copy_nums_and_infos(NeighborConfig {
                max_cycle_size: 10,
                max_flip: 0,
                use_long_cycles: false,
                ignore_cycles_passing_terminal: false,
                use_reducers: false,
            });
            assert_eq!(neighbors.len(), 8);
        }

        {
            let neighbors = dbg.to_neighbor_copy_nums_and_infos(NeighborConfig {
                max_cycle_size: 10,
                max_flip: 2,
                use_long_cycles: false,
                ignore_cycles_passing_terminal: false,
                use_reducers: false,
            });
            assert_eq!(neighbors.len(), 12);
        }
    }
}
