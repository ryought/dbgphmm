//!
//! `SeqGraph` is `DiGraph<N: SeqNode, E: SeqEdge>`
//! seq with copy numbers
//!

use crate::common::CopyNum;
use crate::hmm::params::PHMMParams;
use crate::hmmv2::common::{PEdge, PModel, PNode};
use crate::prob::Prob;
use petgraph::dot::Dot;
use petgraph::graph::DiGraph;
pub use petgraph::graph::{EdgeIndex, NodeIndex};
pub use petgraph::Direction;

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
    fn copy_num(&self) -> Option<CopyNum>;
}

pub struct SeqGraphV2<N: SeqNode, E: SeqEdge>(pub DiGraph<N, E>);

impl<N: SeqNode, E: SeqEdge> SeqGraphV2<N, E> {
    pub fn new(graph: DiGraph<N, E>) -> Self {
        SeqGraphV2(graph)
    }
    /// Get the number of nodes in the SeqGraph
    pub fn node_count(&self) -> usize {
        self.0.node_count()
    }
    /// Get the number of edges in the SeqGraph
    pub fn edge_count(&self) -> usize {
        self.0.edge_count()
    }
}

/*
trait SeqGraph {
    fn total_emittable_copy_num(&self) -> CopyNum;
    fn total_emittable_child_copy_nums(&self, node: NodeIndex) -> CopyNum;
    fn to_phmm_node(&self, node: NodeIndex, total_copy_num: CopyNum) -> PNode;
    fn to_phmm_edge(&self, edge: EdgeIndex) -> PEdge;
    fn to_phmm(&self, param: PHMMParams) -> PModel;
}
*/

impl<N: SeqNode, E: SeqEdge> SeqGraphV2<N, E> {
    /// calculate the sum of copy numbers
    /// of all emittable nodes
    pub fn total_emittable_copy_num(&self) -> CopyNum {
        self.0
            .node_indices()
            .map(|v| {
                let vw = self.0.node_weight(v).unwrap();
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
    pub fn total_emittable_child_copy_nums(&self, node: NodeIndex) -> CopyNum {
        self.0
            .neighbors_directed(node, Direction::Outgoing)
            .map(|child| {
                let child_weight = self.0.node_weight(child).unwrap();
                if child_weight.is_emittable() {
                    child_weight.copy_num()
                } else {
                    0
                }
            })
            .sum()
    }
    /// Convert Node in SimpleSeqGraph into phmm node
    pub fn to_phmm_node(&self, node: NodeIndex, total_copy_num: CopyNum) -> PNode {
        let node_weight = self.0.node_weight(node).unwrap();
        let init_prob =
            Prob::from_prob(node_weight.copy_num() as f64) / Prob::from_prob(total_copy_num as f64);
        PNode::new(
            node_weight.copy_num(),
            init_prob,
            node_weight.is_emittable(),
            node_weight.base(),
        )
    }
    /// Convert Edge in SimpleSeqGraph into phmm edge
    pub fn to_phmm_edge(&self, edge: EdgeIndex) -> PEdge {
        let (parent, child) = self.0.edge_endpoints(edge).unwrap();
        let total_child_copy_num = self.total_emittable_child_copy_nums(parent);
        let child_weight = self.0.node_weight(child).unwrap();
        let trans_prob =
            Prob::from_prob(child_weight.copy_num() as f64 / total_child_copy_num as f64);
        PEdge::new(trans_prob)
    }
    /// convert SimpleSeqGraph to PHMM by ignoreing the edge copy numbers
    pub fn to_phmm(&self, param: PHMMParams) -> PModel {
        let total_copy_num = self.total_emittable_copy_num();
        let graph = self.0.map(
            // node converter
            |v, _| self.to_phmm_node(v, total_copy_num),
            // edge converter
            |e, _| self.to_phmm_edge(e),
        );
        PModel { param, graph }
    }
}

//
// minimum implementations
//

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

impl<N, E> std::fmt::Display for SeqGraphV2<N, E>
where
    N: SeqNode + std::fmt::Display,
    E: SeqEdge + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.0, &[]))
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
}
