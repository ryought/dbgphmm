//!
//! Multi-k edge-centric de Bruijn graph
//!
//! Node: (k-1)-mer
//! Edge: k-mer(s)
//!
//! # Feature
//!
//! * Extend to k+1 edge-centric DBG
//! * Convert to node-centric and PHMM
//! * Not store k-mers in nodes/edges: efficient when k is large
//! * Serializable into GFA sequence graph representation
//! * No generics
//! * Use compact copy number vector
//!
use crate::common::{
    sequence_to_string, CopyNum, Reads, Seq, SeqStyle, Sequence, StyledSequence, NULL_BASE,
};
use crate::dbg::dbg::{Dbg, DbgEdgeBase, DbgNode};
use crate::graph::compact::compact_simple_paths;
use crate::kmer::{
    common::{KmerLike, NullableKmer},
    veckmer::VecKmer,
};
use crate::vector::{DenseStorage, EdgeVec};
use petgraph::graph::{DefaultIx, DiGraph, EdgeIndex, Graph, NodeIndex};
use petgraph::visit::EdgeRef; // for edges_directed
use petgraph::Direction;
use std::convert::From;

///
/// In MultiDbg only edges have copynum.
///
pub type CopyNums = EdgeVec<DenseStorage<CopyNum>>;

///
/// Edge-centric and simple-path-collapsed Dbg structure
/// For large k (k < 10,000)
///
#[derive(Clone, Debug)]
pub struct MultiDbg {
    ///
    /// size of k in de Bruijn graph
    ///
    k: usize,
    ///
    /// Edge-centric de Bruijn graph
    ///
    /// * Edge = k-mer
    ///     emission (last base of k-mer) and corresponding edge in compact
    /// * Node = (k-1)-mer
    ///     terminal or not?
    ///
    full: DiGraph<MultiFullNode, MultiFullEdge>,
    ///
    /// Simple-path-collapsed edge-centric de Bruijn graph
    ///
    /// * Edge = simple path composed of k-mers
    ///     copy_num
    /// * Node = (k-1)-mer
    ///     have same index in graph
    ///
    compact: DiGraph<MultiCompactNode, MultiCompactEdge>,
}

///
/// Edge of MultiDbg full graph
///
#[derive(Clone, Debug)]
pub struct MultiFullEdge {
    ///
    /// emission (last base of k-mer)
    ///
    base: u8,
    ///
    /// corresponding edge index in compact graph
    ///
    edge_in_compact: EdgeIndex,
}

///
/// Node of MultiDbg full graph
///
#[derive(Clone, Debug)]
pub struct MultiFullNode {
    ///
    /// if this node corresponds to k-1mer NNN, then true.
    ///
    is_terminal: bool,
}

///
/// Edge of MultiDbg compact graph
///
#[derive(Clone, Debug)]
pub struct MultiCompactEdge {
    ///
    /// copy number of edge (simple-path or k-mers in full)
    ///
    copy_num: CopyNum,
    ///
    ///
    ///
    edges_in_full: Vec<EdgeIndex>,
}

///
/// Node of MultiDbg compact graph
///
/// Empty struct
///
#[derive(Clone, Debug)]
pub struct MultiCompactNode {}

///
/// Conversion Dbg -> MultiDbg
///
/// Dbg is not simple-path-collapsed, so each node corresponds to a kmer.
///
impl<N: DbgNode, E: DbgEdgeBase> From<Dbg<N, E>> for MultiDbg {
    fn from(dbg: Dbg<N, E>) -> MultiDbg {
        // let graph = dbg.to_edbg_graph(
        //     |km1mer| MultiDbgNode::new(km1mer.is_null()),
        //     |_node, node_weight| {
        //         let seq = vec![node_weight.emission()];
        //         let copy_num = node_weight.copy_num();
        //         MultiDbgEdge::new(seq, copy_num)
        //     },
        // );
        // MultiDbg { k: dbg.k(), graph }
        unimplemented!();
    }
}

//
// Attributes
//

impl MultiDbg {
    ///
    /// Reference of full graph `&DiGraph<MultiFullNode, MultiFullEdge>`
    ///
    pub fn graph_full(&self) -> &DiGraph<MultiFullNode, MultiFullEdge> {
        &self.full
    }
    ///
    /// Reference of compact graph `&DiGraph<(), MultiCompactEdge>`
    ///
    pub fn graph_compact(&self) -> &DiGraph<MultiCompactNode, MultiCompactEdge> {
        &self.compact
    }
    ///
    /// size of k
    ///
    pub fn k(&self) -> usize {
        self.k
    }
    ///
    /// Convert edge in full graph into k-mer
    ///
    pub fn kmer(&self, edge_in_full: EdgeIndex) -> VecKmer {
        unimplemented!();
    }
    ///
    /// Convert edge in compact graph into concated k-mers
    ///
    pub fn kmers(&self, edge_in_compact: EdgeIndex) -> VecKmer {
        unimplemented!();
    }
    ///
    /// Convert node into (k-1)-mer
    /// Concatenate k-1 parental bases
    ///
    pub fn km1mer(&self, node: NodeIndex) -> VecKmer {
        let mut bases = Vec::new();
        let mut node = node;
        let km1 = self.k() - 1;
        while bases.len() < km1 {
            let parent_edge = self
                .graph_full()
                .edges_directed(node, Direction::Incoming)
                .next()
                .unwrap_or_else(|| panic!("no incoming edge"));
            node = parent_edge.source();
            bases.push(parent_edge.weight().base);
        }
        bases.reverse();
        VecKmer::from_bases(&bases)
    }
}

//
// Copy number
//

impl MultiDbg {
    pub fn set_copy_nums(&mut self, copy_nums: &CopyNums) {
        unimplemented!();
    }
    pub fn get_copy_nums(&self) -> CopyNums {
        unimplemented!();
    }
    pub fn neighbor_copy_nums(&self) -> Vec<CopyNums> {
        unimplemented!();
    }
}

//
// Modifying graph structure
//

impl MultiDbg {
    ///
    ///
    ///
    pub fn collapse_all_simple_paths(&self) -> Self {
        let compacted = compact_simple_paths(self.graph_full());
        // let graph = compacted.map(
        //     |_node, weight| weight.clone(),
        //     |_edge, weight| {
        //         MultiDbgEdge::new()
        //     },
        // );
        // MultiDbg { k: self.k(), graph }
        unimplemented!();
    }
    ///
    ///
    ///
    pub fn to_kp1_dbg(self) -> Self {
        unimplemented!();
    }
}

//
// serialize/deserialize
//

impl MultiDbg {
    ///
    /// Dot file with each node/edge shown in Display serialization
    ///
    pub fn to_dot(&self) -> String {
        format!(
            "{}",
            petgraph::dot::Dot::with_config(&self.graph_compact(), &[])
        )
    }
    ///
    ///
    ///
    pub fn to_gfa(&self) -> String {
        unimplemented!();
    }
}

impl std::fmt::Display for MultiCompactEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // todo show edges
        write!(f, "x{}()", self.copy_num)
    }
}

impl std::fmt::Display for MultiFullEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}(E{})",
            self.base as char,
            self.edge_in_compact.index()
        )
    }
}

impl std::fmt::Display for MultiCompactNode {
    fn fmt(&self, _: &mut std::fmt::Formatter) -> std::fmt::Result {
        Ok(())
    }
}

impl std::fmt::Display for MultiFullNode {
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

        for node in multidbg.graph_full().node_indices() {
            let km1mer = multidbg.km1mer(node);
            println!("{:?} {}", node, km1mer);
        }
    }
}
