//!
//! seqgraph
//! seq with copy numbers
//!

use super::common::{PEdge, PGraph, PModel, PNode};
use super::params::PHMMParams;
use crate::common::CopyNum;
use crate::prob::Prob;
use petgraph::dot::Dot;
use petgraph::graph::DiGraph;
pub use petgraph::graph::{EdgeIndex, NodeIndex};
pub use petgraph::Direction;

///
/// SeqGraph is a sequence graph whose node has its own copy number.
///
pub struct SeqGraph<N: SeqNode, E: SeqEdge>(pub DiGraph<N, E>);

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
    /// calculate the sum of copy numbers
    /// of all emittable nodes
    fn total_emittable_copy_num(&self) -> CopyNum {
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
    fn total_emittable_child_copy_nums(&self, node: NodeIndex) -> CopyNum {
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
    /// convert SeqGraph to PHMM by ignoreing the edge copy numbers
    pub fn to_phmm(&self) -> PModel {
        let total_emittable_copy_num = self.total_emittable_copy_num();
        let graph = self.0.map(
            |_, vw| {
                let init_prob = Prob::from_prob(vw.copy_num() as f64)
                    / Prob::from_prob(total_emittable_copy_num as f64);
                PNode::new(vw.copy_num(), init_prob, vw.is_emittable(), vw.base())
            },
            |e, _| {
                let (parent, child) = self.0.edge_endpoints(e).unwrap();
                let total_child_copy_num = self.total_emittable_child_copy_nums(parent);
                let child_w = self.0.node_weight(child).unwrap();
                let trans_prob =
                    Prob::from_prob(child_w.copy_num() as f64 / total_child_copy_num as f64);
                PEdge::new(trans_prob)
            },
        );
        PModel {
            param: PHMMParams::default(),
            graph,
        }
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
        write!(f, "{}", Dot::with_config(&self.0, &[]))
    }
}

//
// mock constructors
//

///
/// Create a linear SeqGraph from a sequence
///
pub fn create_linear_seq_graph(seq: &[u8]) -> SeqGraph<SimpleSeqNode, SimpleSeqEdge> {
    let mut graph = DiGraph::new();
    for i in 0..seq.len() {
        graph.add_node(SimpleSeqNode::new(1, seq[i]));
    }
    for i in 1..seq.len() {
        graph.add_edge(
            NodeIndex::new(i - 1),
            NodeIndex::new(i),
            SimpleSeqEdge::new(None),
        );
    }
    SeqGraph(graph)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_linear_seq_graph_test() {
        let g = create_linear_seq_graph(b"ATCGGCTAGC");
        println!("{}", g);
        let phmm = g.to_phmm();
        println!("{}", phmm);

        for (v, vw) in phmm.nodes() {
            println!("node {:?} {}", v, vw);
            for e in phmm.childs(v) {
                println!("  child {:?}", e);
            }
        }
    }
}
