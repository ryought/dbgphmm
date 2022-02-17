//!
//! `SeqGraph` is `DiGraph<N: SeqNode, E: SeqEdge>`
//! seq with copy numbers
//!

use crate::common::CopyNum;
use crate::hmmv2::common::{PEdge, PModel, PNode};
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;
use petgraph::dot::Dot;
use petgraph::graph::DiGraph;
pub use petgraph::graph::{EdgeIndex, NodeIndex};
pub use petgraph::Direction;

pub trait SeqGraph {
    /// calculate the sum of copy numbers
    /// of all emittable nodes
    fn total_emittable_copy_num(&self) -> CopyNum;
    /// calculate the sum of copy numbers
    /// of all emittable childs of the given node
    fn total_emittable_child_copy_nums(&self, node: NodeIndex) -> CopyNum;
    ///
    /// SeqGraph has consistent copy numbers on nodes?
    ///
    /// * for all nodes, sum of copy numbers of in-edges and out-edges are the same.
    ///
    fn node_copy_nums_is_consistent(&self) -> bool;
    ///
    /// SeqGraph has consistent copy numbers on edges?
    ///
    /// * for all nodes, all of out-edges are either with-copy-numbers or without-copy-numbers
    /// * if copy numbers are set on edges, the sum of copy numbers on edges should be equal to the
    /// copy number of the node.
    ///
    fn edge_copy_nums_is_consistent(&self) -> bool;
    /// Convert Node in SimpleSeqGraph into phmm node
    ///
    fn to_phmm_node(&self, node: NodeIndex, total_copy_num: CopyNum) -> PNode;
    /// Convert Edge in SimpleSeqGraph into phmm edge
    ///
    fn to_phmm_edge(&self, edge: EdgeIndex) -> PEdge;
    /// convert SimpleSeqGraph to PHMM by ignoreing the edge copy numbers
    ///
    fn to_phmm(&self, param: PHMMParams) -> PModel;
}

impl<N: SeqNode, E: SeqEdge> SeqGraph for DiGraph<N, E> {
    /// calculate the sum of copy numbers
    /// of all emittable nodes
    fn total_emittable_copy_num(&self) -> CopyNum {
        self.node_indices()
            .map(|v| {
                let vw = self.node_weight(v).unwrap();
                if vw.is_emittable() {
                    vw.copy_num()
                } else {
                    0
                }
            })
            .sum()
    }
    /// calculate the sum of copy numbers
    /// of all emittable childs of the given node
    fn total_emittable_child_copy_nums(&self, node: NodeIndex) -> CopyNum {
        self.neighbors_directed(node, Direction::Outgoing)
            .map(|child| {
                let child_weight = self.node_weight(child).unwrap();
                if child_weight.is_emittable() {
                    child_weight.copy_num()
                } else {
                    0
                }
            })
            .sum()
    }
    fn node_copy_nums_is_consistent(&self) -> bool {
        // TODO
        true
    }
    fn edge_copy_nums_is_consistent(&self) -> bool {
        self.node_indices().all(|v| {
            let vw = self.node_weight(v).unwrap();
            self.edges_directed(v, Direction::Outgoing)
                .all(|e| e.weight().copy_num().is_some())
        })
    }
    /// Convert Node in SimpleSeqGraph into phmm node
    fn to_phmm_node(&self, node: NodeIndex, total_copy_num: CopyNum) -> PNode {
        let node_weight = self.node_weight(node).unwrap();
        let init_prob = if node_weight.is_emittable() {
            Prob::from_prob(node_weight.copy_num() as f64) / Prob::from_prob(total_copy_num as f64)
        } else {
            Prob::from_prob(0.0)
        };
        PNode::new(
            node_weight.copy_num(),
            init_prob,
            node_weight.is_emittable(),
            node_weight.base(),
        )
    }
    /// Convert Edge in SimpleSeqGraph into phmm edge
    fn to_phmm_edge(&self, edge: EdgeIndex) -> PEdge {
        // TODO if copy_num is assigned?
        let (parent, child) = self.edge_endpoints(edge).unwrap();
        let total_child_copy_num = self.total_emittable_child_copy_nums(parent);
        let child_weight = self.node_weight(child).unwrap();
        let trans_prob = if child_weight.is_emittable() {
            Prob::from_prob(child_weight.copy_num() as f64 / total_child_copy_num as f64)
        } else {
            Prob::from_prob(0.0)
        };
        PEdge::new(trans_prob)
    }
    /// convert SimpleSeqGraph to PHMM by ignoreing the edge copy numbers
    fn to_phmm(&self, param: PHMMParams) -> PModel {
        let total_copy_num = self.total_emittable_copy_num();
        let graph = self.map(
            // node converter
            |v, _| self.to_phmm_node(v, total_copy_num),
            // edge converter
            |e, _| self.to_phmm_edge(e),
        );
        PModel { param, graph }
    }
}

/// a trait that should be satisfyed by nodes in SeqGraph
pub trait SeqNode {
    ///
    /// the copy number of this node
    ///
    fn copy_num(&self) -> CopyNum;
    ///
    /// corresponding base of this node
    ///
    fn base(&self) -> u8;
    ///
    /// This node is emittable or not?
    /// i.e. self.base != 'N'?
    ///
    fn is_emittable(&self) -> bool {
        self.base() != b'N'
    }
}

/// a trait that should be satisfyed by edges in SeqGraph
pub trait SeqEdge {
    ///
    /// the copy number of this edge if assigned
    ///
    fn copy_num(&self) -> Option<CopyNum>;
}

//
// minimum implementations
//

///
/// SimpleSeqGraph is a sequence graph whose node has its own copy number.
///
pub type SimpleSeqGraph = DiGraph<SimpleSeqNode, SimpleSeqEdge>;

pub struct SimpleSeqNode {
    copy_num: CopyNum,
    base: u8,
}

impl SimpleSeqNode {
    pub fn new(copy_num: CopyNum, base: u8) -> SimpleSeqNode {
        SimpleSeqNode { copy_num, base }
    }
}

impl SeqNode for SimpleSeqNode {
    fn copy_num(&self) -> CopyNum {
        self.copy_num
    }
    fn base(&self) -> u8 {
        self.base
    }
}

pub struct SimpleSeqEdge {
    copy_num: Option<CopyNum>,
}

impl SimpleSeqEdge {
    pub fn new(copy_num: Option<CopyNum>) -> SimpleSeqEdge {
        SimpleSeqEdge { copy_num }
    }
}

impl SeqEdge for SimpleSeqEdge {
    fn copy_num(&self) -> Option<CopyNum> {
        self.copy_num
    }
}

//
// Display
//

impl std::fmt::Display for SimpleSeqNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} (x{})", self.base() as char, self.copy_num())
    }
}

impl std::fmt::Display for SimpleSeqEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.copy_num() {
            Some(copy_num) => {
                write!(f, "x{}", copy_num)
            }
            None => {
                write!(f, "")
            }
        }
    }
}

//
// mock constructors
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::mocks::mock_linear;
    #[test]
    fn trait_test() {
        let sg = mock_linear().to_seq_graph();
    }
    #[test]
    fn seq_graph_small() {}
}
