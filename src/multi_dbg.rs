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

use arrayvec::ArrayVec;
use itertools::Itertools;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph::visit::{EdgeRef, IntoNodeReferences};
use petgraph::Direction;
use petgraph_algos::iterators::{ChildEdges, EdgesIterator, NodesIterator, ParentEdges};
use rustflow::min_flow::Flow;

pub mod toy;

///
/// Edge-centric and simple-path-collapsed Dbg structure for large k (k ~ 10,000)
///
#[derive(Clone, Debug)]
pub struct MultiDbg {
    ///
    /// size of k in de Bruijn graph
    ///
    pub k: usize,
    ///
    /// Edge-centric de Bruijn graph
    ///
    /// * Edge = k-mer
    ///     emission (last base of k-mer) and corresponding edge in compact
    /// * Node = (k-1)-mer
    ///     terminal or not?
    ///
    pub full: DiGraph<MultiFullNode, MultiFullEdge>,
    ///
    /// Simple-path-collapsed edge-centric de Bruijn graph
    ///
    /// * Edge = simple path composed of k-mers
    ///     copy_num
    /// * Node = (k-1)-mer
    ///     have same index in graph
    ///
    pub compact: DiGraph<MultiCompactNode, MultiCompactEdge>,
}

///
/// maximum in/out degree (# of parents/childs)
///
pub const MAX_DEGREE: usize = 5;

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
impl<N: DbgNode, E: DbgEdgeBase> std::convert::From<Dbg<N, E>> for MultiDbg {
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

///
/// Graph wrappers and accessor to attributes
///
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
    /// # of edges in compact graph
    ///
    pub fn n_edges_compact(&self) -> usize {
        self.graph_compact().node_count()
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
    ///
    /// Copy number of edge in full graph
    ///
    pub fn copy_num(&self, edge_in_full: EdgeIndex) -> CopyNum {
        self.graph_full()[edge_in_full].copy_num
    }
    ///
    /// Base (emission) of edge in full graph
    ///
    pub fn base(&self, edge_in_full: EdgeIndex) -> u8 {
        self.graph_full()[edge_in_full].base
    }
}

///
/// Bridge functions between full and compact
///
impl MultiDbg {
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
    /// Determine copy number of edge in compact graph
    ///
    /// = copy number of one of the corresponding edges in full graph
    ///
    fn copy_num_of_edge_in_compact(&self, edge_in_compact: EdgeIndex) -> CopyNum {
        unimplemented!();
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

///
/// CopyNums
///
pub type CopyNums = Flow<CopyNum>;

///
/// # Copy number related functions
///
/// `CopyNums` (=`Flow<CopyNum>`) stores copy numbers of all edges in compact.
///
///
impl MultiDbg {
    ///
    /// For all nodes
    ///
    pub fn is_copy_nums_valid(&self) -> bool {
        unimplemented!();
    }
    pub fn genome_size(&self) -> CopyNum {
        unimplemented!();
    }
    pub fn set_copy_nums(&mut self, copy_nums: &CopyNums) {
        assert_eq!(copy_nums.len(), self.n_edges_compact());
        for (e, _, _, _) in self.edges_compact() {
            // let c = copy_nums[e];
        }
        unimplemented!();
    }
    pub fn get_copy_nums(&self) -> CopyNums {
        let mut copy_nums = CopyNums::new(self.n_edges_compact(), 0);
        for (edge, _, _, edge_weight) in self.edges_compact() {
            // TODO
            // copy_nums[edge] =
            // let c = copy_nums[e];
        }
        copy_nums
    }
    pub fn neighbor_copy_nums(&self) -> Vec<CopyNums> {
        unimplemented!();
    }
    ///
    ///
    ///
    pub fn guess_copy_num_of_kp1_edge(
        &self,
        node: NodeIndex,
        edge_in: EdgeIndex,
        edge_out: EdgeIndex,
    ) -> CopyNum {
        let copy_num_in = self.copy_num(edge_in);
        let copy_num_outs: ArrayVec<CopyNum, MAX_DEGREE> =
            self.childs_full(node).map(|(_, _, w)| w.copy_num).collect();
        let out_index = self
            .childs_full(node)
            .position(|(edge, _, _)| edge == edge_out)
            .expect("edge_out is not a child of node");
        Self::guess_copy_num(copy_num_in, &copy_num_outs, out_index)
    }
    /// Assign copy number to the pair of edges incoming/outgoing from the same node.
    /// (used in `guess_copy_num_of_kp1_edge`.)
    ///
    /// Give the same amount of copy number to all >0x outgoing edge in order of copy numbers of
    /// outgoing edge.
    ///
    /// ```
    /// use dbgphmm::multi_dbg::MultiDbg;
    /// assert_eq!(MultiDbg::guess_copy_num(10, &[1, 1, 0], 0), 5);
    /// assert_eq!(MultiDbg::guess_copy_num(10, &[1, 1, 0], 1), 5);
    /// assert_eq!(MultiDbg::guess_copy_num(10, &[1, 1, 0], 2), 0);
    ///
    /// assert_eq!(MultiDbg::guess_copy_num(9, &[1, 1, 0], 0), 5);
    /// assert_eq!(MultiDbg::guess_copy_num(9, &[1, 1, 0], 1), 4);
    /// assert_eq!(MultiDbg::guess_copy_num(9, &[1, 1, 0], 2), 0);
    ///
    /// assert_eq!(MultiDbg::guess_copy_num(1, &[0, 1, 1], 0), 0);
    /// assert_eq!(MultiDbg::guess_copy_num(1, &[0, 1, 2], 1), 0);
    /// assert_eq!(MultiDbg::guess_copy_num(1, &[0, 1, 2], 2), 1);
    ///
    /// assert_eq!(MultiDbg::guess_copy_num(5, &[2, 0, 9, 7], 0), 1);
    /// assert_eq!(MultiDbg::guess_copy_num(5, &[2, 0, 9, 7], 1), 0);
    /// assert_eq!(MultiDbg::guess_copy_num(5, &[2, 0, 9, 7], 2), 2);
    /// assert_eq!(MultiDbg::guess_copy_num(5, &[2, 0, 9, 7], 3), 2);
    ///
    /// assert_eq!(MultiDbg::guess_copy_num(5, &[1], 0), 5);
    /// ```
    pub fn guess_copy_num(
        copy_num_in: CopyNum,
        copy_num_outs: &[CopyNum],
        out_index: usize,
    ) -> CopyNum {
        assert!(out_index < copy_num_outs.len());
        let n_outs = copy_num_outs.iter().filter(|&c| *c > 0).count();

        // Determine out_index is k-th element in copy_num_outs?
        let mut v: ArrayVec<(usize, CopyNum), MAX_DEGREE> =
            copy_num_outs.iter().copied().enumerate().collect();
        // reverse sort by copy_num
        v.sort_by(|(_, copy_num_a), (_, copy_num_b)| copy_num_b.cmp(copy_num_a));
        let rank = v
            .iter()
            .position(|(index, _copy_num)| *index == out_index)
            .unwrap();

        let copy_num_out = copy_num_outs[out_index];
        if copy_num_out == 0 {
            0
        } else {
            let mut copy_num = copy_num_in / n_outs;
            if rank < copy_num_in % n_outs {
                copy_num += 1;
            }
            copy_num
        }
    }
}

///
/// Modifying graph structure
///
impl MultiDbg {
    ///
    /// Construct compact graph from full graph by collapsing
    ///
    pub fn construct_compact_from_full(
        full: &DiGraph<MultiFullNode, MultiFullEdge>,
    ) -> DiGraph<MultiCompactNode, MultiCompactEdge> {
        compact_simple_paths_for_targeted_nodes(full, |node_weight| !node_weight.is_terminal).map(
            |_node, _| MultiCompactNode::new(),
            |_edge, edge_weight| {
                let edges = edge_weight.into_iter().map(|(edge, _)| *edge).collect();
                MultiCompactEdge::new(edges)
            },
        )
    }
    /// Extend k to k+1.
    ///
    ///
    pub fn to_kp1_dbg(&self) -> Self {
        // create full
        let full = self.to_node_centric_graph(
            |_, _| MultiFullNode::new(false), // to_node
            || MultiFullNode::new(true),      // to_terminal_node
            // to_edge
            |e_in, e_out, node| {
                // base: e_out's base
                let base = self.base(e_out);
                // copy_num: determined by `guess_copy_num`. distribute e_in's copy_num
                let copy_num = self.guess_copy_num_of_kp1_edge(node, e_in, e_out);
                MultiFullEdge::new(base, copy_num)
            },
            // to_terminal_edge
            |e| {
                let w = &self.graph_full()[e];
                MultiFullEdge::new(w.base, w.copy_num)
            },
        );

        // create compact from full
        let compact = Self::construct_compact_from_full(&full);

        MultiDbg {
            k: self.k() + 1,
            full,
            compact,
        }
    }
    /// Construct node centric (full) graph G'
    ///
    /// ```text
    /// node v in G' == edge (k-mer) v in G (id will be conserved)
    /// edge from node v to node w in G' == (target of edge v is source of edge w in G)
    /// ```
    ///
    /// Use in k+1 extension and PHMM conversion.
    ///
    /// to_node: `(edge_index, edge_weight)`
    /// to_edge: `(edge_index_source, edge_index_target, center_node_index)`
    ///
    /// # Example
    ///
    /// k=4
    /// ```text
    /// +---+ ATCC
    /// |ATC+--+
    /// +---+  +-->+---+     +---+
    ///            |TCC+---->|CCA|
    /// +---+  +-->+---+     +---+
    /// |GTC+--+        TCCA
    /// +---+ GTCC
    /// ```
    ///
    /// k=5
    /// ```text
    ///        ATCCA
    /// +----+
    /// |ATCC+---+
    /// +----+   +-->+----+
    ///              |TCCA|
    /// +----+   +-->+----+
    /// |GTCC+---+
    /// +----+
    ///        GTCCA
    /// ```
    ///
    /// # Procedure
    ///
    /// ## (1) Add nodes in G' corresponding to each edge in G
    ///
    /// ```text
    /// G:
    ///        e
    /// --> v ---> w -->
    ///
    /// G':
    /// -> to_node(e) ->
    /// ```
    ///
    /// ## (2) For each node in G, add an edge between a pair of incoming edge v and outgoing edge w.
    ///
    /// ```text
    /// G:
    ///  e1         e2
    /// ----> node ---->
    ///
    /// G':
    ///       to_edge(e1,e2,node)
    ///  e1 ----------------------> e2
    /// ```
    ///
    /// If the node is terminal, add new node to represent terminal.
    ///
    /// ```text
    /// G:
    ///       GNNN        NNNA
    /// GNN  -----> NNN  -----> NNA
    ///
    /// G':
    /// -> GNNN --> NNNN --> NNNA ->
    /// ```
    ///
    pub fn to_node_centric_graph<N, E, FN, FTN, FE, FTE>(
        &self,
        to_node: FN,
        to_terminal_node: FTN,
        to_edge: FE,
        to_terminal_edge: FTE,
    ) -> DiGraph<N, E>
    where
        FN: Fn(EdgeIndex, &MultiFullEdge) -> N,
        FTN: Fn() -> N,
        FE: Fn(EdgeIndex, EdgeIndex, NodeIndex) -> E,
        FTE: Fn(EdgeIndex) -> E,
    {
        let mut graph = DiGraph::new();

        // convert edge in G into node in G'
        //
        let to_node_index = |edge: EdgeIndex| NodeIndex::new(edge.index());

        // (1) add a node corresponding to each edge
        //
        for (edge, _, _, edge_weight) in self.edges_full() {
            let node = graph.add_node(to_node(edge, edge_weight));
            assert_eq!(node.index(), edge.index(), "edge index is corrupted");
        }

        // (2) for each node
        for (node, node_weight) in self.nodes_full() {
            if node_weight.is_terminal {
                let terminal_node = graph.add_node(to_terminal_node());
                for (e, _, _) in self.parents_full(node) {
                    let v = to_node_index(e);
                    graph.add_edge(v, terminal_node, to_terminal_edge(e));
                }
                for (e, _, _) in self.childs_full(node) {
                    let v = to_node_index(e);
                    graph.add_edge(terminal_node, v, to_terminal_edge(e));
                }
            } else {
                for (e1, _, _) in self.parents_full(node) {
                    let v1 = to_node_index(e1);
                    for (e2, _, _) in self.childs_full(node) {
                        let v2 = to_node_index(e2);
                        graph.add_edge(v1, v2, to_edge(e1, e2, node));
                    }
                }
            }
        }

        graph
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
        write!(f, "{}({}x)", self.base as char, self.copy_num)
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
