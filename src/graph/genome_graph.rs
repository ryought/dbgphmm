//!
//! `GenomeGraph`
//! Node: linear sequence
//! Edge: its adjacency
//!
use super::seq_graph::SeqGraph;
use crate::common::{CopyNum, Sequence};
use petgraph::dot::Dot;
use petgraph::graph::DiGraph;

/// GenomeGraph
pub struct GenomeGraph(DiGraph<GenomeNode, GenomeEdge>);

/// Node in GenomeGraph is a copy-number-assigned sequence fragment
pub struct GenomeNode {
    /// sequence fragment of this node
    seq: Sequence,
    /// copy number of the sequence
    copy_num: CopyNum,
}

impl GenomeNode {
    fn new(seq: &[u8], copy_num: CopyNum) -> Self {
        GenomeNode {
            seq: seq.to_vec(),
            copy_num,
        }
    }
}

pub struct GenomeEdge {
    /// copy number of the transition
    /// it assumes random transition if `None`
    copy_num: Option<CopyNum>,
}

impl GenomeEdge {
    fn new(copy_num: Option<CopyNum>) -> Self {
        GenomeEdge { copy_num }
    }
}

impl std::fmt::Display for GenomeNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let seq = std::str::from_utf8(&self.seq).unwrap();
        write!(f, "{} (x{})", seq, self.copy_num)
    }
}

impl std::fmt::Display for GenomeEdge {
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

impl std::fmt::Display for GenomeGraph {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.0, &[]))
    }
}

/*
impl GenomeGraph {
    fn to_seq_graph(&self) -> SeqGraph {
        unimplemented!();
    }
}
*/

pub fn mock() -> GenomeGraph {
    let mut g = DiGraph::new();
    let v1 = g.add_node(GenomeNode::new(b"ATTCGATCGTCG", 1));
    let v2 = g.add_node(GenomeNode::new(b"ATCGATG", 1));
    let v3 = g.add_node(GenomeNode::new(b"TTCGAT", 2));
    g.add_edge(v1, v3, GenomeEdge::new(None));
    g.add_edge(v2, v3, GenomeEdge::new(None));
    GenomeGraph(g)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genome_graph_mock1() {
        let g = mock();
        println!("{}", g);
    }
}
