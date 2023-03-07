//!
//! Multi-k edge-centric de Bruijn graph
//!
//! Node: (k-1)-mer
//! Edge: k-mer(s)
//!
//! # Feature
//!
//! * Easy to k+1 extension (edge-centric)
//! * Less number of nodes/edges (collapses simple path)
//! * Memory-efficient when k is large
//! * Convertable to node-centric when converting to PHMM
//! * Serializable into GFA sequence graph representation
//! * No generics
//!
//! # Todos
//! * From/Into Dbg and PHMM
//! * calculate score
//!
use crate::common::{CopyNum, Reads, Seq, SeqStyle, Sequence, StyledSequence, NULL_BASE};
use crate::kmer::veckmer::VecKmer;
use crate::segment::Segment;
use petgraph::graph::{DefaultIx, DiGraph, EdgeIndex, Graph, NodeIndex};

struct MultiDbg {
    ///
    /// k
    ///
    k: usize,
    ///
    /// backend graph of Dbg
    ///
    graph: DiGraph<MultiDbgNode, MultiDbgEdge>,
}

struct MultiDbgEdge {
    segment: Segment,
    copy_num: CopyNum,
}

struct MultiDbgNode {
    km1mer: Segment,
    is_terminal: bool,
}

impl MultiDbg {
    ///
    ///
    ///
    pub fn to_kp1_dbg(self) -> Self {
        unimplemented!();
    }
    pub fn set_copy_nums(&mut self) {}
    pub fn get_copy_nums(&self) {}
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f() {}
}
