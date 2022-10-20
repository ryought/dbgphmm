//!
//! copy number enumeration with cycle basis
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase, NodeCopyNums};
use super::edge_centric::impls::{SimpleEDbgEdge, SimpleEDbgNode};
use super::edge_centric::EDbgNode;
use crate::graph::cycle_space::CycleSpace;
use crate::graph::spanning_tree::spanning_tree;
use crate::kmer::kmer::{Kmer, KmerLike};
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
}
impl<K: KmerLike> Iterator for CopyNumsIterator<K> {
    type Item = NodeCopyNums;
    fn next(&mut self) -> Option<Self::Item> {
        // update copy_nums using each cycle (for both +1/-1 directions) if resulting copy_nums
        // vector is valid.
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
}
