//!
//! copy number enumeration with cycle basis
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase, EdgeCopyNums, NodeCopyNums};
use super::edge_centric::compact::{
    compacted_flow_into_original_flow, into_compacted_flow, SimpleCompactedEDbgEdge,
};
use super::edge_centric::impls::{SimpleEDbgEdge, SimpleEDbgNode};
use super::edge_centric::{EDbgEdge, EDbgEdgeBase, EDbgEdgeMin, EDbgNode};
use crate::dbg::hashdbg_v2::HashDbg;
use crate::graph::cycle::{
    apply_cycle_with_dir, to_cycle_with_dir, Cycle, CycleWithDir, SimpleCycle,
};
use crate::graph::cycle_space::CycleSpace;
use crate::graph::spanning_tree::spanning_tree;
use crate::hist::{DiscreteDistribution, Hist};
use crate::kmer::common::concat_overlapping_kmers;
use crate::kmer::kmer::{Kmer, KmerLike};
use crate::min_flow::base::FlowEdgeBase;
use crate::min_flow::enumerate_neighboring_flows;
use crate::min_flow::residue::{
    total_changes, update_info_to_cycle_with_dir, ResidueDirection, UpdateInfo,
};
use crate::prob::Prob;
use crate::utils::all_same_value;
use crate::vector::{DenseStorage, NodeVec};
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

#[derive(Clone, Debug)]
pub struct CopyNumsUpdateInfo<K: KmerLike> {
    /// segment info (cyclic path in residue graph as kmer and direction segments)
    segments: Vec<(K, ResidueDirection)>,
    ///
    genome_size_change: i32,
    /// size of cycle (= total number of edges in cycle) in compacted/uncompacted graph
    cycle_size: usize,
    /// in uncompacted graph it will be equal to `cycle_size`.
    n_kmer_changed: usize,
}

impl<K: KmerLike> std::fmt::Display for CopyNumsUpdateInfo<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}(dG={},L={},n={})",
            self.segments
                .iter()
                .map(|(kmer, dir)| format!("{}{}", kmer, dir))
                .join(","),
            self.genome_size_change,
            self.cycle_size,
            self.n_kmer_changed,
        )
    }
}

///
/// - k: k-size of dbg
///
pub fn edges_to_kmer<N, E: EDbgEdgeMin>(
    edbg: &DiGraph<N, E>,
    edges: &[EdgeIndex],
    k: usize,
) -> E::Kmer {
    let kmers: Vec<_> = edges
        .iter()
        .map(|&edge| edbg.edge_weight(edge).unwrap().kmer().clone())
        .collect();
    concat_overlapping_kmers(kmers, k - 1)
}

impl<K: KmerLike> CopyNumsUpdateInfo<K> {
    pub fn empty() -> Self {
        CopyNumsUpdateInfo {
            segments: vec![],
            genome_size_change: 0,
            cycle_size: 0,
            n_kmer_changed: 0,
        }
    }
    pub fn from_uncompacted(
        k: usize,
        graph: &DiGraph<SimpleEDbgNode<K>, SimpleEDbgEdge<K>>,
        update_info: &UpdateInfo,
    ) -> Self {
        let cycle_with_dir = update_info_to_cycle_with_dir(update_info);
        let segments: Vec<_> = cycle_with_dir
            .collapse_dir()
            .into_iter()
            .map(|(mut edges, is_rev)| {
                if is_rev {
                    edges.reverse();
                    (edges_to_kmer(graph, &edges, k), ResidueDirection::Down)
                } else {
                    (edges_to_kmer(graph, &edges, k), ResidueDirection::Up)
                }
            })
            .collect();
        let genome_size_change = total_changes(update_info.iter().map(|(_, dir)| *dir));
        CopyNumsUpdateInfo {
            segments,
            genome_size_change,
            cycle_size: update_info.len(),
            n_kmer_changed: update_info.len(),
        }
    }
    pub fn from_compacted(
        k: usize,
        graph: &DiGraph<SimpleEDbgNode<K>, SimpleCompactedEDbgEdge<K>>,
        update_info: &UpdateInfo,
    ) -> Self {
        let cycle_with_dir = update_info_to_cycle_with_dir(update_info);
        let segments: Vec<_> = cycle_with_dir
            .collapse_dir()
            .into_iter()
            .map(|(mut edges, is_rev)| {
                if is_rev {
                    edges.reverse();
                    (edges_to_kmer(graph, &edges, k), ResidueDirection::Down)
                } else {
                    (edges_to_kmer(graph, &edges, k), ResidueDirection::Up)
                }
            })
            .collect();
        let genome_size_change = update_info
            .iter()
            .map(|&(edge, dir)| {
                let n_edges_in_uncompacted = graph.edge_weight(edge).unwrap().origin_edges().len();
                dir.int() * n_edges_in_uncompacted as i32
            })
            .sum();
        let n_kmer_changed = update_info
            .iter()
            .map(|&(edge, _)| {
                let n_edges_in_uncompacted = graph.edge_weight(edge).unwrap().origin_edges().len();
                n_edges_in_uncompacted
            })
            .sum();
        CopyNumsUpdateInfo {
            segments,
            genome_size_change,
            n_kmer_changed,
            cycle_size: update_info.len(),
        }
    }
}

fn have_zero_one_change(a: &EdgeCopyNums, b: &EdgeCopyNums) -> bool {
    assert_eq!(a.len(), b.len());
    for i in 0..(a.len()) {
        let v = EdgeIndex::new(i);
        if (a[v] == 0 && b[v] > 0) || (a[v] > 0 && b[v] == 0) {
            return true;
        }
    }
    false
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
    /// wrapper of `neighbor_copy_nums_fast`
    pub fn neighbor_copy_nums_fast(&self) -> Vec<NodeCopyNums> {
        self.neighbor_copy_nums_fast_with_info()
            .into_iter()
            .map(|(copy_nums, _)| copy_nums)
            .collect()
    }
    /// wrapper of `neighbor_copy_nums_fast_compact`
    pub fn neighbor_copy_nums_fast_compact(&self, max_depth: usize) -> Vec<NodeCopyNums> {
        self.neighbor_copy_nums_fast_compact_with_info(max_depth, false)
            .into_iter()
            .map(|(copy_nums, _)| copy_nums)
            .collect()
    }
    ///
    /// use Johnson1975
    ///
    pub fn neighbor_copy_nums_fast_with_info(
        &self,
    ) -> Vec<(NodeCopyNums, CopyNumsUpdateInfo<N::Kmer>)> {
        // convert to edbg, residue graph
        let edbg = self.to_edbg();
        let network = edbg.graph.map(
            |_, _| (),
            |_, weight| {
                let copy_num = weight.copy_num();
                FlowEdgeBase::new(copy_num.saturating_sub(1), copy_num.saturating_add(1), 0.0)
            },
        );
        let copy_num = self.to_node_copy_nums().switch_index();
        // enumerate all cycles
        enumerate_neighboring_flows(&network, &copy_num, None)
            .into_iter()
            .map(|(flow, update_info)| {
                (
                    flow.switch_index(),
                    CopyNumsUpdateInfo::from_uncompacted(self.k(), &edbg.graph, &update_info),
                )
            })
            .collect()
    }
    ///
    /// use Johnson1975 on Compacted Edbg
    ///
    /// * If ignore_high_copys is true, only new copy numbers that have 1x->0x or 0x->1x change is
    /// returned.
    ///
    pub fn neighbor_copy_nums_fast_compact_with_info(
        &self,
        max_depth: usize,
        ignore_high_copys: bool,
    ) -> Vec<(NodeCopyNums, CopyNumsUpdateInfo<N::Kmer>)> {
        let graph = self.to_compact_edbg_graph();
        // println!("{}", petgraph::dot::Dot::with_config(&graph, &[]));
        let network = graph.map(
            |_, _| (),
            |_, weight| {
                let copy_num = weight.copy_num();
                FlowEdgeBase::new(copy_num.saturating_sub(1), copy_num.saturating_add(1), 0.0)
            },
        );
        let copy_num = into_compacted_flow(&graph, &self.to_node_copy_nums().switch_index());
        // enumerate all cycles
        enumerate_neighboring_flows(&network, &copy_num, Some(max_depth))
            .into_iter()
            .filter(|(new_copy_num, _)| {
                if ignore_high_copys {
                    have_zero_one_change(&copy_num, new_copy_num)
                } else {
                    true
                }
            })
            .map(|(flow, update_info)| {
                let flow_in_original =
                    compacted_flow_into_original_flow(self.n_nodes(), &graph, &flow);
                (
                    flow_in_original.switch_index(),
                    CopyNumsUpdateInfo::from_compacted(self.k(), &graph, &update_info),
                )
            })
            .collect()
    }
    ///
    ///
    pub fn to_kmer_distribution(
        &self,
        neighbors: &[(NodeCopyNums, Prob)],
    ) -> Vec<DiscreteDistribution> {
        self.nodes()
            .map(|(node, weight)| {
                let copy_nums: Vec<_> = neighbors.iter().map(|(cn, p)| cn[node]).collect();
                let hist = Hist::from(&copy_nums);
                let copy_nums_with_prob: Vec<_> =
                    neighbors.iter().map(|(cn, p)| (cn[node], *p)).collect();
                DiscreteDistribution::from_occurs(&copy_nums_with_prob)
            })
            .collect()
    }
    ///
    /// Remove all nodes `v` such that `P(c_v=0) >= p_0` and current copy_num is 0.
    ///
    pub fn purge_zero_copy_with_high_prob_kmer(
        &mut self,
        dds: &[DiscreteDistribution],
        p_0: Prob,
    ) -> usize {
        assert_eq!(dds.len(), self.n_nodes());
        let n_before = self.graph.node_count();
        self.graph.retain_nodes(|g, v| {
            let is_likely_error_kmer = dds[v.index()].p_x(0) >= p_0;
            let is_zero_copy = g.node_weight(v).unwrap().copy_num() == 0;
            !(is_likely_error_kmer && is_zero_copy)
        });
        let n_after = self.graph.node_count();
        n_before - n_after
    }
    ///
    ///
    pub fn to_copy_num_expected_vector(
        &self,
        dds: &[DiscreteDistribution],
    ) -> NodeVec<DenseStorage<f64>> {
        assert_eq!(dds.len(), self.n_nodes());
        let mut ret = NodeVec::new(self.n_nodes(), 0.0);
        for (node, _) in self.nodes() {
            ret[node] = dds[node.index()].mean();
        }
        ret
    }
    ///
    /// check the variance of copy_num of each kmer
    ///
    pub fn inspect_kmer_variance(
        &self,
        neighbors: &[(NodeCopyNums, Prob)],
        copy_nums_true: &NodeCopyNums,
        read_count: &HashDbg<N::Kmer>,
    ) {
        self.inspect_kmer_variance_with_comment(neighbors, copy_nums_true, read_count, |_| {
            String::new()
        })
    }
    ///
    /// check the variance of copy_num of each kmer
    ///
    pub fn inspect_kmer_variance_with_comment<F>(
        &self,
        neighbors: &[(NodeCopyNums, Prob)],
        copy_nums_true: &NodeCopyNums,
        read_count: &HashDbg<N::Kmer>,
        comment: F,
    ) where
        F: Fn(NodeIndex) -> String,
    {
        let k = self.k();
        let print_header = || {
            println!(
                "#K k={}\tkmer\tnode_id\tcurrent_copy_num\ttrue_copy_num\tread_count\tp(copy_num=copy_num_true)\tp(copy_num=0)\tprobs\tdegree_info(in,out)\thist\tcopy_nums\tcomment",
                k
            );
        };
        print_header();
        let kmer_distributions = self.to_kmer_distribution(neighbors);
        for (node, weight) in self
            .nodes()
            .sorted_by_key(|&(node, _)| copy_nums_true[node])
        {
            let copy_nums: Vec<_> = neighbors.iter().map(|(cn, p)| cn[node]).collect();
            let hist = Hist::from(&copy_nums);
            // let copy_nums_with_prob: Vec<_> = neighbors
            //     .iter()
            //     .map(|(cn, p)| (cn[node], p.to_value()))
            //     .collect();
            let copy_num_true = copy_nums_true[node];
            println!(
                "K\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}",
                k,
                weight.kmer(),
                node.index(),
                weight.copy_num(),
                copy_num_true,
                read_count.get(weight.kmer()),
                kmer_distributions[node.index()]
                    .p_x(copy_num_true)
                    .to_value(),
                kmer_distributions[node.index()].p_x(0).to_value(),
                kmer_distributions[node.index()],
                format!("({},{})", self.in_degree(node), self.out_degree(node)),
                hist,
                // copy_nums_with_prob,
                comment(node),
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
        assert_eq!(copy_nums_v2.len(), copy_nums.len());

        println!("c3");
        let copy_nums_v3 = dbg.neighbor_copy_nums_fast_compact(100);
        for copy_num in copy_nums_v3.iter() {
            println!("c3={}", copy_num);
        }
        assert!(is_equal_as_set(&copy_nums, &copy_nums_v3));
        assert_eq!(copy_nums_v3.len(), copy_nums.len());
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
    #[test]
    fn neighbor_copy_nums_update_info_test() {
        let dbg = mock_intersection_small();

        // uncompacted
        let neighbors = dbg.neighbor_copy_nums_fast_with_info();
        for (copy_num, info) in neighbors.iter() {
            println!("c={} info={}", copy_num, info);
        }

        // compacted
        let neighbors = dbg.neighbor_copy_nums_fast_compact_with_info(100, false);
        for (copy_num, info) in neighbors.iter() {
            println!("c={} info={}", copy_num, info);
        }
    }
}
