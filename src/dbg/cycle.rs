//!
//! copy number enumeration with cycle basis
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase, EdgeCopyNums, NodeCopyNums};
use super::edge_centric::impls::{SimpleEDbgEdge, SimpleEDbgNode};
use super::edge_centric::EDbgNode;
use crate::graph::cycle::{
    apply_cycle_with_dir, to_cycle_with_dir, Cycle, CycleWithDir, SimpleCycle,
};
use crate::graph::cycle_space::CycleSpace;
use crate::graph::spanning_tree::spanning_tree;
use crate::hist::{get_normalized_probs, Hist};
use crate::kmer::kmer::{Kmer, KmerLike};
use crate::prob::Prob;
use itertools::Itertools;
use petgraph::dot::Dot;
use petgraph::graph::{NodeIndex, UnGraph};

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

///
/// * undirected edge-centric dbg
/// * CycleSpace of the graph
///
pub struct CopyNumsIterator<K: KmerLike> {
    edbg: UndirectedEdbg<K>,
    space: CycleSpace,
    // queue: Vec<NodeCopyNums>,
}
impl<K: KmerLike> Iterator for CopyNumsIterator<K> {
    type Item = NodeCopyNums;
    fn next(&mut self) -> Option<Self::Item> {
        // (1) pop a candidate update cycle from cycle space
        // (2) update copy_nums using each cycle (for both +1/-1 directions)
        // if resulting copy_nums vector is valid.
        unimplemented!();
    }
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
    pub fn iter_neighbor_copy_nums(&self) -> CopyNumsIterator<N::Kmer> {
        // (1) create undirected edge-centric-dbg.
        let edbg = self.to_undirected_edbg_graph();
        // (2) enumerate cycles in the undirected graph using CycleSpace iterator.
        let null_node = get_null_node(&edbg);
        let st = spanning_tree(&edbg, null_node);
        let basis = st.cycle_basis_list(&edbg);
        let space = CycleSpace::new(basis);
        CopyNumsIterator { edbg, space }
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
                        // println!("c+={}", c);
                        ret.push((c, cycle_increase.clone()));
                    } else {
                        // println!("c+=no");
                    }

                    // -1 along cycle
                    let cycle_decrease = cycle_increase.reverse_dir();
                    let decreased = apply_cycle_with_dir(&copy_nums, &cycle_decrease);
                    if decreased.is_some() {
                        let c = decreased.unwrap().switch_index();
                        // println!("c-={}", c);
                        ret.push((c, cycle_decrease));
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
                .map(|(x, p)| format!("p(x={})={}", x, p.to_value()))
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
        // let segments = cycle.collapse_dir();
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
        format!("len={}", cycle.edges().len())
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::super::mocks::*;
    use super::*;
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
    }
}
