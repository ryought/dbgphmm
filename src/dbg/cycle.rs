//!
//! copy number enumeration with cycle basis
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase, EdgeCopyNums, NodeCopyNums};
use super::edge_centric::impls::{SimpleEDbgEdge, SimpleEDbgNode};
use super::edge_centric::{EDbgEdge, EDbgEdgeBase, EDbgNode};
use crate::graph::cycle::{
    apply_cycle_with_dir, to_cycle_with_dir, Cycle, CycleWithDir, SimpleCycle,
};
use crate::graph::cycle_space::CycleSpace;
use crate::graph::spanning_tree::spanning_tree;
use crate::hist::{get_normalized_probs, Hist};
use crate::kmer::kmer::{Kmer, KmerLike};
use crate::min_flow::residue::enumerate_neighboring_flows_in_residue;
use crate::prob::Prob;
use crate::utils::all_same_value;
use fnv::FnvHashSet as HashSet;
use itertools::Itertools;
use petgraph::dot::Dot;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex, UnGraph};

pub type UndirectedEdbg<K> = UnGraph<SimpleEDbgNode<K>, SimpleEDbgEdge<K>>;

///
/// get NodeIndex of null-km1mer (NNN).
///
fn get_null_node<K: KmerLike>(edbg: &UndirectedEdbg<K>) -> NodeIndex {
    edbg.node_indices()
        .find(|&node| {
            let weight = edbg.node_weight(node).unwrap();
            weight.km1mer().is_null()
        })
        .expect("edbg does not contain null km1mer")
}

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    ///
    ///
    fn to_undirected_edbg_graph(&self) -> UndirectedEdbg<N::Kmer> {
        self.to_edbg_graph(
            |km1mer| SimpleEDbgNode::new(km1mer.clone()),
            |node, weight| {
                let kmer = weight.kmer().clone();
                let copy_num = weight.copy_num();
                SimpleEDbgEdge::new(kmer, copy_num, node)
            },
        )
    }
    ///
    /// Get an iterator of copy_nums (NodeVec of CopyNum) neighboring the current copy_nums of
    /// self.
    ///
    /// ```text
    /// for copy_nums in dbg.iter_neighbor_copy_nums() {
    /// }
    /// ```
    ///
    pub fn iter_neighbor_copy_nums(&self) {
        unimplemented!();
    }
    pub fn neighbor_copy_nums(&self) -> Vec<NodeCopyNums> {
        self.neighbor_copy_nums_and_cycles()
            .into_iter()
            .map(|(ncn, _)| ncn)
            .collect()
    }
    ///
    /// List up all neighboring copy numbers
    /// **Heavy function**
    ///
    pub fn neighbor_copy_nums_and_cycles(&self) -> Vec<(NodeCopyNums, CycleWithDir)> {
        // (1) create undirected edge-centric-dbg.
        let edbg = self.to_undirected_edbg_graph();
        // (2) enumerate cycles in the undirected graph using CycleSpace iterator.
        let null_node = get_null_node(&edbg);
        let st = spanning_tree(&edbg, null_node);
        let basis = st.cycle_basis_list(&edbg);
        let space = CycleSpace::new(basis);

        let mut ret = Vec::new();
        let mut is_stored: HashSet<NodeCopyNums> = HashSet::default();
        let edbg_directed = self.to_edbg();
        let copy_nums = self.to_node_copy_nums().switch_index();
        // println!("{}", Dot::with_config(&edbg_directed.graph, &[]));
        for simple_cycle in space {
            // println!("simple_cycle={}", simple_cycle);
            match simple_cycle.to_cycle(&edbg) {
                Some(cycle) => {
                    // +1 along cycle
                    let cycle_increase = to_cycle_with_dir(&edbg_directed.graph, &cycle);
                    let increased = apply_cycle_with_dir(&copy_nums, &cycle_increase);
                    if increased.is_some() {
                        let c = increased.unwrap().switch_index();
                        if !is_stored.contains(&c) {
                            // c is new copy_nums vector
                            // println!("c+={}", c);
                            ret.push((c.clone(), cycle_increase.clone()));
                            is_stored.insert(c);
                        }
                    } else {
                        // println!("c+=no");
                    }

                    // -1 along cycle
                    let cycle_decrease = cycle_increase.reverse_dir();
                    let decreased = apply_cycle_with_dir(&copy_nums, &cycle_decrease);
                    if decreased.is_some() {
                        let c = decreased.unwrap().switch_index();
                        if !is_stored.contains(&c) {
                            // println!("c-={}", c);
                            ret.push((c.clone(), cycle_decrease));
                            is_stored.insert(c);
                        }
                    } else {
                        // println!("c-=no");
                    }
                }
                None => {}
            };
        }

        ret
    }
    ///
    /// use Johnson1975
    ///
    pub fn neighbor_copy_nums_fast(&self) -> Vec<NodeCopyNums> {
        // convert to edbg, residue graph
        let rg = self.to_residue_edbg();
        let copy_num = self.to_node_copy_nums().switch_index();
        // enumerate all cycles
        enumerate_neighboring_flows_in_residue(&rg, &copy_num, None)
            .into_iter()
            .map(|(flow, _)| flow.switch_index())
            .collect()
    }
    ///
    /// use Johnson1975 on Compacted Edbg
    ///
    pub fn neighbor_copy_nums_fast_compact(&self) -> Vec<NodeCopyNums> {
        let graph = self.to_compact_edbg_graph();

        // convert to edbg, residue graph
        let rg = self.to_residue_edbg();
        let copy_num = self.to_node_copy_nums().switch_index();
        // enumerate all cycles
        enumerate_neighboring_flows_in_residue(&rg, &copy_num, None)
            .into_iter()
            .map(|(flow, _)| flow.switch_index())
            .collect()
    }
    ///
    /// check the variance of copy_num of each kmer
    ///
    pub fn inspect_kmer_variance(&self, neighbors: &[(NodeCopyNums, Prob)]) {
        let print_header = || {
            println!("#K kmer\tnode_id\tcurrent_copy_num\tprobs\tcopy_nums");
        };

        print_header();
        for (node, weight) in self.nodes() {
            let copy_nums: Vec<_> = neighbors.iter().map(|(cn, p)| cn[node]).collect();
            let hist = Hist::from(&copy_nums);

            let copy_nums_with_prob: Vec<_> =
                neighbors.iter().map(|(cn, p)| (cn[node], *p)).collect();
            let normalized = get_normalized_probs(&copy_nums_with_prob);
            let txt_normalized = normalized
                .iter()
                .map(|(x, p)| format!("p(x={})={:.5}", x, p.to_value()))
                .join(",");
            println!(
                "K\t{}\t{}\t{}\t{}\t{}\t{:?}",
                weight.kmer(),
                node.index(),
                weight.copy_num(),
                hist,
                txt_normalized,
                copy_nums,
            );
        }
        print_header();
    }
    pub fn summarize_cycle_with_dir(&self, cycle: &CycleWithDir) -> String {
        let segments = cycle.collapse_dir();
        // segments
        //     .iter()
        //     .map(|(edges, is_rev)| {
        //         format!(
        //             "{}{}",
        //             edges
        //                 .iter()
        //                 .map(|edge| self.kmer(NodeIndex::new(edge.index())))
        //                 .join(","),
        //             is_rev
        //         )
        //     })
        //     .join(",")
        // format!("len={},edges={}", cycle.edges().len(), cycle)
        format!(
            "len={},segments={}",
            cycle.edges().len(),
            segments
                .iter()
                .map(|(edges, is_rev)| if *is_rev {
                    format!("{}-", edges.len())
                } else {
                    format!("{}+", edges.len())
                })
                .join(",")
        )
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::super::mocks::*;
    use super::*;
    use crate::utils::is_equal_as_set;
    use petgraph::dot::Dot;
    #[test]
    fn basic_dbg_to_undirected() {
        // let dbg = mock_simple();
        let dbg = mock_intersection();
        let ug = dbg.to_undirected_edbg_graph();
        println!("{}", Dot::with_config(&ug, &[]));

        let null_node = get_null_node(&ug);
        println!("{:?}", null_node);
        // assert_eq!(null_node, NodeIndex::new(7));

        let st = spanning_tree(&ug, null_node);
        let basis = st.cycle_basis_list(&ug);
        for (i, b) in basis.iter().enumerate() {
            println!("#{} {}", i, b);
        }
    }
    ///
    /// abbr of `NodeCopyNums::from_slice(xs, 0)` used in the test `neighbor_copy_nums_test`
    ///
    fn to_nc(xs: &[usize]) -> NodeCopyNums {
        NodeCopyNums::from_slice(xs, 0)
    }
    #[test]
    fn neighbor_copy_nums_test() {
        let dbg = mock_intersection_small();
        println!("c0={}", dbg.to_node_copy_nums());
        let copy_nums = dbg.neighbor_copy_nums();

        // debug
        for copy_num in copy_nums.iter() {
            println!("c={}", copy_num);
        }

        assert_eq!(
            copy_nums,
            vec![
                to_nc(&[1, 2, 1, 1, 1, 0, 0, 0, 2, 1, 1, 1, 2, 1, 2, 1, 0, 1]),
                to_nc(&[1, 0, 1, 1, 1, 2, 2, 2, 0, 1, 1, 1, 0, 1, 0, 1, 2, 1]),
                to_nc(&[1, 2, 2, 2, 1, 1, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 1, 1]),
                to_nc(&[1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1]),
                to_nc(&[2, 2, 1, 1, 2, 1, 1, 1, 2, 2, 1, 1, 2, 2, 2, 1, 1, 2]),
                to_nc(&[0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0]),
                to_nc(&[1, 1, 2, 2, 1, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1]),
                to_nc(&[1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1]),
                to_nc(&[2, 1, 1, 1, 2, 2, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 2, 2]),
                to_nc(&[0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0]),
                to_nc(&[2, 1, 0, 0, 2, 1, 1, 1, 1, 2, 0, 0, 1, 2, 1, 0, 1, 2]),
                to_nc(&[0, 1, 2, 2, 0, 1, 1, 1, 1, 0, 2, 2, 1, 0, 1, 2, 1, 0]),
            ]
        );

        let copy_nums_v2 = dbg.neighbor_copy_nums_fast();
        for copy_num in copy_nums_v2.iter() {
            println!("c2={}", copy_num);
        }

        assert!(is_equal_as_set(&copy_nums, &copy_nums_v2));
    }
    #[test]
    fn compact_edbg_01() {
        {
            let dbg = mock_simple();
            println!("{}", Dot::with_config(&dbg.to_edbg().graph, &[]));
            println!("{}", Dot::with_config(&dbg.to_compact_edbg_graph(), &[]));
        }

        {
            let dbg = mock_intersection();
            println!("{}", Dot::with_config(&dbg.to_edbg().graph, &[]));
            println!("{}", Dot::with_config(&dbg.to_compact_edbg_graph(), &[]));
        }
    }
}
