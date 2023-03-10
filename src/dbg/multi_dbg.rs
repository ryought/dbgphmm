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
use crate::graph::compact::compact_simple_paths_for_targeted_nodes;
use crate::kmer::{
    common::{KmerLike, NullableKmer},
    veckmer::VecKmer,
};
use crate::vector::{DenseStorage, EdgeVec};
use petgraph::graph::{DefaultIx, DiGraph, EdgeIndex, Graph, NodeIndex};
use petgraph::visit::EdgeRef; // for edges_directed
use petgraph::Direction;
use std::convert::From;

use itertools::Itertools;
use petgraph_algos::iterators::{ChildEdges, EdgesIterator, NodesIterator, ParentEdges};

pub mod toy;

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
/// Corresponds to k-mer
///
#[derive(Clone, Debug)]
pub struct MultiFullEdge {
    ///
    /// emission (last base of k-mer)
    ///
    base: u8,
    ///
    /// copy number of k-mer
    ///
    copy_num: CopyNum,
    // ///
    // /// corresponding edge index in compact graph
    // ///
    // edge_in_compact: EdgeIndex,
}

impl MultiFullEdge {
    pub fn new(base: u8, copy_num: CopyNum) -> MultiFullEdge {
        MultiFullEdge { base, copy_num }
    }
}

///
/// Node of MultiDbg full graph
///
/// Corresponds to (k-1)-mer
///
#[derive(Clone, Debug)]
pub struct MultiFullNode {
    ///
    /// if this node corresponds to k-1mer NNN, then true.
    ///
    is_terminal: bool,
}

impl MultiFullNode {
    pub fn new(is_terminal: bool) -> MultiFullNode {
        MultiFullNode { is_terminal }
    }
}

///
/// Edge of MultiDbg compact graph
///
#[derive(Clone, Debug)]
pub struct MultiCompactEdge {
    // ///
    // /// copy number of edge (simple-path or k-mers in full)
    // ///
    // copy_num: CopyNum,
    ///
    /// corresponding edges in full graph
    ///
    /// Ordering of edges is source to terminal.
    /// For example, if below four edges in full is collapsed into an edge in compact, then `[e0,e1,e2,e3]`.
    ///
    /// ```text
    /// e0 e1 e2 e3
    /// -->-->-->-->
    /// ```
    ///
    edges_in_full: Vec<EdgeIndex>,
}

impl MultiCompactEdge {
    pub fn new(edges_in_full: Vec<EdgeIndex>) -> Self {
        Self { edges_in_full }
    }
}

///
/// Node of MultiDbg compact graph
///
/// Empty struct
///
#[derive(Clone, Debug)]
pub struct MultiCompactNode {}

impl MultiCompactNode {
    pub fn new() -> Self {
        Self {}
    }
}

///
/// Conversion Dbg (w/ copy number) -> MultiDbg
///
impl<N: DbgNode, E: DbgEdgeBase> From<Dbg<N, E>> for MultiDbg {
    fn from(dbg: Dbg<N, E>) -> MultiDbg {
        let full = dbg.to_edbg_graph(
            |km1mer| MultiFullNode::new(km1mer.is_null()),
            |_node, node_weight| MultiFullEdge::new(node_weight.emission(), node_weight.copy_num()),
        );

        // create compact from
        let compact = Self::construct_compact_from_full(&full);

        MultiDbg {
            k: dbg.k(),
            full,
            compact,
        }
    }
}

//
// Graph wrappers
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
    /// Iterator of all nodes in the graph
    /// Item is (node: NodeIndex, node_weight: &MultiFullNode)
    pub fn nodes_full(&self) -> NodesIterator<MultiFullNode> {
        NodesIterator::new(self.graph_full())
    }
    /// Iterator of all edges in the graph
    /// Item is (edge: EdgeIndex, source: NodeIndex, target: NodeIndex, edge_weight: &MultiFullEdge)
    pub fn edges_full(&self) -> EdgesIterator<MultiFullEdge> {
        EdgesIterator::new(self.graph_full())
    }
    /// Iterator of childs of the node
    /// Item is (edge: EdgeIndex, child: NodeIndex, edge_weight: &MultiFullEdge)
    pub fn childs_full(&self, node: NodeIndex) -> ChildEdges<MultiFullEdge> {
        ChildEdges::new(self.graph_full(), node)
    }
    /// Iterator of parents of the node
    /// Item is (edge: EdgeIndex, parent: NodeIndex, edge_weight: &MultiFullEdge)
    pub fn parents_full(&self, node: NodeIndex) -> ParentEdges<MultiFullEdge> {
        ParentEdges::new(self.graph_full(), node)
    }
    /// Iterator of all nodes in the graph
    /// Item is (node: NodeIndex, node_weight: &MultiCompactNode)
    pub fn nodes_compact(&self) -> NodesIterator<MultiCompactNode> {
        NodesIterator::new(self.graph_compact())
    }
    /// Iterator of all edges in the graph
    /// Item is (edge: EdgeIndex, source: NodeIndex, target: NodeIndex, edge_weight: &MultiCompactEdge)
    pub fn edges_compact(&self) -> EdgesIterator<MultiCompactEdge> {
        EdgesIterator::new(self.graph_compact())
    }
    /// Iterator of childs of the node
    /// Item is (edge: EdgeIndex, child: NodeIndex, edge_weight: &MultiCompactEdge)
    pub fn childs_compact(&self, node: NodeIndex) -> ChildEdges<MultiCompactEdge> {
        ChildEdges::new(self.graph_compact(), node)
    }
    /// Iterator of parents of the node
    /// Item is (edge: EdgeIndex, parent: NodeIndex, edge_weight: &MultiCompactEdge)
    pub fn parents_compact(&self, node: NodeIndex) -> ParentEdges<MultiCompactEdge> {
        ParentEdges::new(self.graph_compact(), node)
    }
}

///
/// Determine kmer of each node/edge
///
impl MultiDbg {
    ///
    /// Convert edge in full graph into k-mer
    ///
    pub fn kmer_full(&self, edge_in_full: EdgeIndex) -> VecKmer {
        let (source, _) = self.graph_full().edge_endpoints(edge_in_full).unwrap();
        let base = self.graph_full().edge_weight(edge_in_full).unwrap().base;
        let km1mer = self.km1mer_full(source);
        km1mer.extend_last(base)
    }
    ///
    /// Convert node in full graph into (k-1)-mer
    /// (Concatenate k-1 parental bases)
    ///
    pub fn km1mer_full(&self, node_in_full: NodeIndex) -> VecKmer {
        let mut bases = Vec::new();
        let mut node = node_in_full;
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
    ///
    /// Determine corresponding node in full
    ///
    fn map_node_in_compact_to_full(&self, node_in_compact: NodeIndex) -> NodeIndex {
        //
        // pick a outgoing edge from the node.
        // map the edge in compact into a simple path in full.
        // the source node of first edge in the simple path is the corresponding node in full.
        let (_, _, ew) = self.childs_compact(node_in_compact).next().unwrap();
        let first_child_edge = ew.edges_in_full[0];
        let (node_in_full, _) = self.graph_full().edge_endpoints(first_child_edge).unwrap();
        node_in_full
    }
    ///
    /// Convert node in compact graph into (k-1)-mer
    ///
    pub fn km1mer_compact(&self, node_in_compact: NodeIndex) -> VecKmer {
        let node_in_full = self.map_node_in_compact_to_full(node_in_compact);
        self.km1mer_full(node_in_full)
    }
    ///
    /// Convert edge in compact graph into concated k-mers
    ///
    /// (kmer of edge in compact) = (km1mer of source node) + (bases in edge)
    ///
    pub fn kmer_compact(&self, edge_in_compact: EdgeIndex) -> VecKmer {
        let (source_node, _) = self
            .graph_compact()
            .edge_endpoints(edge_in_compact)
            .unwrap();
        let mut kmer = self.km1mer_compact(source_node);

        for &edge_in_full in self
            .graph_compact()
            .edge_weight(edge_in_compact)
            .unwrap()
            .edges_in_full
            .iter()
        {
            let base = self.graph_full().edge_weight(edge_in_full).unwrap().base;
            kmer = kmer.into_extend_last(base);
        }

        kmer
    }
}

//
// Copy number
//

impl MultiDbg {
    pub fn is_copy_nums_valid(&self) -> bool {
        unimplemented!();
    }
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
    /// Construct compact graph from full graph by collapsing
    ///
    pub fn construct_compact_from_full(
        full: &DiGraph<MultiFullNode, MultiFullEdge>,
    ) -> DiGraph<MultiCompactNode, MultiCompactEdge> {
        compact_simple_paths_for_targeted_nodes(full, |node_weight| !node_weight.is_terminal).map(
            |node, node_weight| MultiCompactNode::new(),
            |edge, edge_weight| {
                let edges = edge_weight.into_iter().map(|(edge, _)| *edge).collect();
                MultiCompactEdge::new(edges)
            },
        )
    }
    ///
    ///
    ///
    pub fn collapse_all_simple_paths(&self) -> Self {
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
            "Full:\n{}Compact:\n{}",
            petgraph::dot::Dot::with_config(&self.graph_full(), &[]),
            petgraph::dot::Dot::with_config(&self.graph_compact(), &[]),
        )
    }
    ///
    /// Debug output each node/edge with kmer
    ///
    pub fn show_graph_with_kmer(&self) {
        println!("Full:");
        for (node, weight) in self.nodes_full() {
            println!("v{}\t{}\t{}", node.index(), weight, self.km1mer_full(node));
        }
        for (edge, s, t, weight) in self.edges_full() {
            println!(
                "e{}\tv{}\tv{}\t{}\t{}",
                edge.index(),
                s.index(),
                t.index(),
                weight,
                self.kmer_full(edge)
            );
        }

        println!("Compact:");
        for (node, weight) in self.nodes_compact() {
            println!(
                "v{}\t{}\t{}",
                node.index(),
                weight,
                self.km1mer_compact(node)
            );
        }
        for (edge, s, t, weight) in self.edges_compact() {
            println!(
                "e{}\tv{}\tv{}\t{}\t{}",
                edge.index(),
                s.index(),
                t.index(),
                weight,
                self.kmer_compact(edge)
            );
        }
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
        write!(
            f,
            "{}",
            self.edges_in_full
                .iter()
                .map(|e| format!("e{}", e.index()))
                .join(",")
        )
    }
}

impl std::fmt::Display for MultiFullEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}(x{})", self.base as char, self.copy_num)
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
            let km1mer = multidbg.km1mer_full(node);
            println!("{:?} {}", node, km1mer);
        }
    }
}
