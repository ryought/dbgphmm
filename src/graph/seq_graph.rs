//!
//! `SeqGraph`
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

///
/// SeqGraph is a sequence graph whose node has its own copy number.
///
pub struct SeqGraph<N: SeqNode, E: SeqEdge> {
    graph: DiGraph<N, E>,
    total_emittable_copy_num: CopyNum,
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
    fn copy_num(&self) -> Option<CopyNum>;
}

impl<N: SeqNode, E: SeqEdge> SeqGraph<N, E> {
    /// Create `SeqGraph` from `DiGraph<N: SeqNode, E: SeqEdge>`
    /// with precalculation of total_emittable_copy_num
    pub fn new(graph: DiGraph<N, E>) -> Self {
        let total_emittable_copy_num = SeqGraph::total_emittable_copy_num(&graph);
        SeqGraph {
            graph,
            total_emittable_copy_num,
        }
    }
    /// Get the number of nodes in the SeqGraph
    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }
    /// Get the number of edges in the SeqGraph
    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }
}

//
// Conversion between SeqGraph and GenomeGraph
//
impl<N: SeqNode, E: SeqEdge> SeqGraph<N, E> {
    /// calculate the sum of copy numbers
    /// of all emittable nodes
    fn total_emittable_copy_num(graph: &DiGraph<N, E>) -> CopyNum {
        graph
            .node_indices()
            .map(|v| {
                let vw = graph.node_weight(v).unwrap();
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
    fn total_emittable_child_copy_nums(graph: &DiGraph<N, E>, node: NodeIndex) -> CopyNum {
        graph
            .neighbors_directed(node, Direction::Outgoing)
            .map(|child| {
                let child_weight = graph.node_weight(child).unwrap();
                if child_weight.is_emittable() {
                    child_weight.copy_num()
                } else {
                    0
                }
            })
            .sum()
    }
    /// Convert Node in Seqgraph into phmm node
    fn to_phmm_node(&self, node: NodeIndex) -> PNode {
        let node_weight = self.graph.node_weight(node).unwrap();
        let init_prob = Prob::from_prob(node_weight.copy_num() as f64)
            / Prob::from_prob(self.total_emittable_copy_num as f64);
        PNode::new(
            node_weight.copy_num(),
            init_prob,
            node_weight.is_emittable(),
            node_weight.base(),
        )
    }
    /// Convert Edge in Seqgraph into phmm edge
    fn to_phmm_edge(&self, edge: EdgeIndex) -> PEdge {
        let (parent, child) = self.graph.edge_endpoints(edge).unwrap();
        let total_child_copy_num = SeqGraph::total_emittable_child_copy_nums(&self.graph, parent);
        let child_weight = self.graph.node_weight(child).unwrap();
        let trans_prob =
            Prob::from_prob(child_weight.copy_num() as f64 / total_child_copy_num as f64);
        PEdge::new(trans_prob)
    }
    /// convert SeqGraph to PHMM by ignoreing the edge copy numbers
    pub fn to_phmm(&self, param: PHMMParams) -> PModel {
        let graph = self.graph.map(
            // node converter
            |v, _| self.to_phmm_node(v),
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

impl<N, E> std::fmt::Display for SeqGraph<N, E>
where
    N: SeqNode + std::fmt::Display,
    E: SeqEdge + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.graph, &[]))
    }
}

//
// mock constructors
//

#[cfg(test)]
mod tests {
    use super::*;
}
