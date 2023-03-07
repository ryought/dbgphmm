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
use crate::common::{
    sequence_to_string, CopyNum, Reads, Seq, SeqStyle, Sequence, StyledSequence, NULL_BASE,
};
use crate::dbg::dbg::{Dbg, DbgEdgeBase, DbgNode};
use crate::kmer::{
    common::{KmerLike, NullableKmer},
    veckmer::VecKmer,
};
use petgraph::graph::{DefaultIx, DiGraph, EdgeIndex, Graph, NodeIndex};
use std::convert::From;

///
/// Edge-centric and simple-path-collapsed Dbg structure
///
/// # Todos
///
/// ## Constructor
/// * convert from Dbg
///
/// ## Attributes
/// * convert edge/node into kmer/kp1mer
///
/// ## Copy number
/// * neighbor
///
/// ## Modifying graph structure
/// * extend into k+1 dbg
/// * collapse simple path
///
/// ## Dumping
/// * serialize/deserialize: GFA and string representation
///
#[derive(Clone, Debug)]
pub struct MultiDbg {
    ///
    /// k
    ///
    k: usize,
    ///
    /// Backend graph of Dbg using petgraph
    /// Node is MultiDbgNode and edge is MultiDbgEdge
    ///
    graph: DiGraph<MultiDbgNode, MultiDbgEdge>,
}

///
/// Edge of MultiDbg
///
#[derive(Clone, Debug)]
pub struct MultiDbgEdge {
    seq: Sequence,
    copy_num: CopyNum,
}

impl MultiDbgEdge {
    ///
    /// Constructor
    ///
    pub fn new(seq: Sequence, copy_num: CopyNum) -> Self {
        MultiDbgEdge { seq, copy_num }
    }
    ///
    /// Reference to copy_num
    ///
    pub fn copy_num(&self) -> CopyNum {
        self.copy_num
    }
    ///
    /// Reference to seq
    ///
    pub fn seq(&self) -> &Sequence {
        &self.seq
    }
}

///
/// Node of MultiDbg
///
#[derive(Clone, Debug)]
pub struct MultiDbgNode {
    ///
    /// if this node corresponds to k-1mer NNN, then true.
    ///
    is_terminal: bool,
}

impl MultiDbgNode {
    ///
    /// Constructor
    ///
    pub fn new(is_terminal: bool) -> Self {
        MultiDbgNode { is_terminal }
    }
    ///
    /// Reference to is_terminal
    ///
    pub fn is_terminal(&self) -> bool {
        self.is_terminal
    }
}

///
/// Conversion Dbg -> MultiDbg
///
/// Dbg is not simple-path-collapsed, so each node corresponds to a kmer.
///
impl<N: DbgNode, E: DbgEdgeBase> From<Dbg<N, E>> for MultiDbg {
    fn from(dbg: Dbg<N, E>) -> MultiDbg {
        let graph = dbg.to_edbg_graph(
            |km1mer| MultiDbgNode::new(km1mer.is_null()),
            |_node, node_weight| {
                let seq = vec![node_weight.emission()];
                let copy_num = node_weight.copy_num();
                MultiDbgEdge::new(seq, copy_num)
            },
        );
        MultiDbg { k: dbg.k(), graph }
    }
}

impl MultiDbg {
    ///
    /// Reference to backend graph (petgraph::DiGraph<MultiDbgNode, MultiDbgEdge>)
    ///
    pub fn graph(&self) -> &DiGraph<MultiDbgNode, MultiDbgEdge> {
        &self.graph
    }
    ///
    /// Convert edge into concatenated kmers
    ///
    pub fn kmers(&self, edge: EdgeIndex) -> VecKmer {
        unimplemented!();
    }
    ///
    /// Convert node into (k-1)mer
    ///
    pub fn km1mer(&self, node: NodeIndex) -> VecKmer {
        unimplemented!();
    }
    ///
    ///
    ///
    pub fn to_kp1_dbg(self) -> Self {
        unimplemented!();
    }
    pub fn set_copy_nums(&mut self) {}
    pub fn get_copy_nums(&self) {}
    pub fn neighbor_copy_nums(&self) {}
}

//
// serialize/deserialize
//

impl MultiDbg {
    ///
    /// Dot file with each node/edge shown in Display serialization
    ///
    pub fn to_dot(&self) -> String {
        format!("{}", petgraph::dot::Dot::with_config(&self.graph, &[]))
    }
    ///
    ///
    ///
    pub fn to_gfa(&self) -> String {
        unimplemented!();
    }
}

// impl std::fmt::Display for MultiDbg {
//     fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
//         unimplemented!();
//     }
// }

impl std::fmt::Display for MultiDbgEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}(x{})", sequence_to_string(&self.seq), self.copy_num)
    }
}

impl std::fmt::Display for MultiDbgNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "is_terminal={}", self.is_terminal)
    }
}

///
/// Mock example MultiDbg definitions for debug
///
mod mocks {
    use super::*;
    ///
    /// Circular `ATCTCCG` in k=4
    ///
    /// `ATCT`
    /// `TCTC`
    /// `CTCC`
    /// `TCCG`
    /// `CCGA`
    /// `CGAT`
    /// `GATC`
    ///
    pub fn multidbg_circular() -> MultiDbg {
        unimplemented!();
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks as dbgmocks;

    #[test]
    fn convert_from_dbg() {
        let dbg = dbgmocks::mock_intersection_small();
        println!("{}", dbg);
        let multidbg: MultiDbg = dbg.into();
        println!("{}", multidbg.to_dot());
    }
}
