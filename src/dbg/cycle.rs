//!
//! copy number enumeration with cycle basis
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase, EdgeCopyNums, NodeCopyNums};
use super::edge_centric::impls::{SimpleEDbgEdge, SimpleEDbgNode};
use super::edge_centric::EDbgNode;
use crate::graph::cycle::{apply_cycle, Cycle, SimpleCycle};
use crate::graph::cycle_space::CycleSpace;
use crate::graph::spanning_tree::spanning_tree;
use crate::kmer::kmer::{Kmer, KmerLike};
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
    ///
    /// List up all neighboring copy numbers
    /// **Heavy function**
    ///
    pub fn neighbor_copy_nums(&self) -> Vec<NodeCopyNums> {
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
                    let increased = apply_cycle(&edbg_directed.graph, &copy_nums, &cycle, false);
                    if increased.is_some() {
                        let c = increased.unwrap().switch_index();
                        // println!("c+={}", c);
                        ret.push(c);
                    } else {
                        // println!("c+=no");
                    }
                    // -1 along cycle
                    let decreased = apply_cycle(&edbg_directed.graph, &copy_nums, &cycle, true);
                    if decreased.is_some() {
                        let c = decreased.unwrap().switch_index();
                        // println!("c-={}", c);
                        ret.push(c);
                    } else {
                        // println!("c-=no");
                    }
                }
                None => {}
            };
        }

        ret
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
    #[test]
    fn neighbor_test() {
        let dbg = mock_intersection_small();
        println!("c0={}", dbg.to_node_copy_nums());
        for copy_num in dbg.neighbor_copy_nums() {
            println!("c={}", copy_num);
        }
    }
}
