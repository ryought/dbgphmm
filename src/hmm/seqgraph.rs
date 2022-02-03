//!
//! seqgraph
//! seq with copy numbers
//!

use crate::common::CopyNum;
use petgraph::dot::Dot;
use petgraph::graph::DiGraph;
pub use petgraph::graph::{EdgeIndex, NodeIndex};

pub struct SeqGraph(pub DiGraph<SeqNode, SeqEdge>);

pub struct SeqNode {
    copy_num: CopyNum,
    base: u8,
}

impl SeqNode {
    pub fn new(copy_num: CopyNum, base: u8) -> SeqNode {
        SeqNode { copy_num, base }
    }
}

pub struct SeqEdge {
    copy_num: Option<CopyNum>,
}

impl SeqEdge {
    pub fn new(copy_num: Option<CopyNum>) -> SeqEdge {
        SeqEdge { copy_num }
    }
}

//
// Display
//

impl std::fmt::Display for SeqNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} (x{})", self.base as char, self.copy_num)
    }
}

impl std::fmt::Display for SeqEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.copy_num {
            Some(copy_num) => {
                write!(f, "x{}", copy_num)
            }
            None => {
                write!(f, "")
            }
        }
    }
}

impl std::fmt::Display for SeqGraph {
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
pub fn create_linear_seq_graph(seq: &[u8]) -> SeqGraph {
    let mut graph = DiGraph::new();
    for i in 0..seq.len() {
        graph.add_node(SeqNode::new(1, seq[i]));
    }
    for i in 1..seq.len() {
        graph.add_edge(NodeIndex::new(i - 1), NodeIndex::new(i), SeqEdge::new(None));
    }
    SeqGraph(graph)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_linear_seq_graph_test() {
        let g = create_linear_seq_graph(b"ATCGGCTAGCT");
        println!("{}", g);
    }
}
