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
use crate::dbg::hashdbg_v2::HashDbg;
use crate::graph::compact::compact_simple_paths_for_targeted_nodes;
use crate::graph::utils::degree_stats;
use crate::kmer::{
    common::{KmerLike, NullableKmer},
    veckmer::VecKmer,
};

use arrayvec::ArrayVec;
use fnv::FnvHashMap as HashMap;
use itertools::Itertools;
use petgraph::graph::{DefaultIx, DiGraph, EdgeIndex, NodeIndex};
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
    ///
    /// Base (emission) is null (`b'n'`)?
    ///
    pub fn is_null_base(&self) -> bool {
        self.base == NULL_BASE
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
    pub fn edges_in_full(&self) -> &[EdgeIndex] {
        &self.edges_in_full
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
/// Conversion HashDbg -> MultiDbg
///
impl<K: KmerLike> std::convert::From<HashDbg<K>> for MultiDbg {
    fn from(hashdbg: HashDbg<K>) -> MultiDbg {
        unimplemented!();
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
    /// # of nodes in full graph
    ///
    pub fn n_nodes_full(&self) -> usize {
        self.graph_full().node_count()
    }
    ///
    /// # of nodes in compact graph
    ///
    pub fn n_nodes_compact(&self) -> usize {
        self.graph_compact().node_count()
    }
    ///
    /// # of edges in full graph
    ///
    pub fn n_edges_full(&self) -> usize {
        self.graph_full().edge_count()
    }
    ///
    /// # of edges in compact graph
    ///
    pub fn n_edges_compact(&self) -> usize {
        self.graph_compact().edge_count()
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
    /// Copy number of node in full graph = sum of incoming/outgoing edges of the node.
    ///
    pub fn copy_num_of_node(&self, node_in_full: NodeIndex) -> CopyNum {
        self.childs_full(node_in_full)
            .map(|(_, _, edge_weight)| edge_weight.copy_num)
            .sum()
    }
    ///
    /// Base (emission) of edge in full graph
    ///
    pub fn base(&self, edge_in_full: EdgeIndex) -> u8 {
        self.graph_full()[edge_in_full].base
    }
    ///
    /// Get terminal node (= (k-1)-mer NNN) in full
    ///
    pub fn terminal_node_full(&self) -> Option<NodeIndex> {
        self.nodes_full()
            .find(|(_, node_weight)| node_weight.is_terminal)
            .map(|(node, _)| node)
    }
    ///
    /// count the number of nodes with (in_degree, out_degree) in graph full
    ///
    pub fn degree_stats(&self) -> HashMap<(usize, usize), usize> {
        degree_stats(self.graph_full())
    }
    ///
    /// Count # of ambiguous nodes whose in_deg > 1 and out_deg > 1.
    ///
    pub fn n_ambiguous_node(&self) -> usize {
        let mut n = 0;
        for node in self.graph_full().node_indices() {
            let in_degree = self.parents_full(node).count();
            let out_degree = self.childs_full(node).count();
            if in_degree > 1 && out_degree > 1 {
                n += 1;
            }
        }
        n
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
    /// Convert a path (= sequence of edges) in compact into a path in full.
    ///
    pub fn path_full_from_path_compact(&self, path_compact: &[EdgeIndex]) -> Vec<EdgeIndex> {
        let mut path_full = Vec::new();
        for &edge_in_compact in path_compact {
            for &edge_in_full in self.graph_compact()[edge_in_compact].edges_in_full() {
                path_full.push(edge_in_full);
            }
        }
        path_full
    }
}

///
/// Path and styled seq conversion related
///
impl MultiDbg {
    ///
    ///
    ///
    pub fn path_to_styled_seqs(&self, path_full: &[EdgeIndex]) -> Vec<StyledSequence> {
        let mut seqs = Vec::new();
        let mut state: Option<(Vec<u8>, NodeIndex, SeqStyle)> = None;
        let terminal_node = self.terminal_node_full();

        for &e in path_full {
            let (s, t) = self
                .graph_full()
                .edge_endpoints(e)
                .expect("edge in path is not in graph");

            // beginning of new seq
            if state.is_none() {
                let style = if terminal_node.is_some() && s == terminal_node.unwrap() {
                    SeqStyle::Linear
                } else {
                    SeqStyle::Circular
                };
                state = Some((Vec::new(), s, style));
            }

            let (mut seq, first_node, style) = state.unwrap();
            let weight = &self.graph_full()[e];
            if weight.is_null_base() {
                assert!(!style.is_circular());
            } else {
                seq.push(weight.base);
            }

            if t == first_node {
                // end of this seq
                seqs.push(StyledSequence::new(seq, style));
                state = None
            } else {
                // iterate again
                state = Some((seq, first_node, style));
            }
        }
        seqs
    }
    ///
    /// Generate Euler circuit path of full graph that traverses all edges
    ///
    pub fn get_euler_circuit(&self) -> Vec<EdgeIndex> {
        let mut path = Vec::new();

        // start from terminal node if exists
        let start_node = self.terminal_node_full().unwrap_or(NodeIndex::new(0));

        let mut node = start_node;
        let mut copy_nums_remain = self.get_copy_nums_full();

        loop {
            match self
                .childs_full(node)
                .find(|(edge, _, _)| copy_nums_remain[*edge] > 0)
            {
                Some((edge, child, _)) => {
                    path.push(edge);
                    copy_nums_remain[edge] -= 1;
                    node = child;
                }
                None => {
                    if node == start_node {
                        break;
                    } else {
                        panic!("euler traverse not found");
                    }
                }
            }
        }

        path
    }
    ///
    /// used in `get_euler_circuit` to create copy_nums_remain
    ///
    /// CopyNums of edges in full is abnormal
    ///
    fn get_copy_nums_full(&self) -> CopyNums {
        let mut copy_nums = CopyNums::new(self.n_edges_full(), 0);
        for (edge, _, _, edge_weight) in self.edges_full() {
            copy_nums[edge] = edge_weight.copy_num;
        }
        copy_nums
    }
    ///
    ///
    pub fn to_styled_seqs(&self) -> Vec<StyledSequence> {
        self.path_to_styled_seqs(&self.get_euler_circuit())
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
    ///
    /// Convert edge in compact graph into sequence (emission bases)
    ///
    pub fn seq_compact(&self, edge_in_compact: EdgeIndex) -> Vec<u8> {
        let mut seq = Vec::new();

        for &edge_in_full in self
            .graph_compact()
            .edge_weight(edge_in_compact)
            .unwrap()
            .edges_in_full()
            .iter()
        {
            let base = self.graph_full().edge_weight(edge_in_full).unwrap().base;
            seq.push(base);
        }

        seq
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
    /// Check if copy numbers assigned to all edges in full is valid (satisfying flow-consistency).
    ///
    /// For all nodes v, (sum of copy_num of incoming edges) == (sum of copy_num of outgoing edges)
    ///
    pub fn is_copy_nums_valid(&self) -> bool {
        self.nodes_full().all(|(node, _)| {
            let copy_num_in: CopyNum = self.parents_full(node).map(|(_, _, w)| w.copy_num).sum();
            let copy_num_out: CopyNum = self.childs_full(node).map(|(_, _, w)| w.copy_num).sum();
            copy_num_in == copy_num_out
        })
    }
    ///
    /// Genome size is a sum of copy number of all edges with non-null emission.
    ///
    pub fn genome_size(&self) -> CopyNum {
        self.edges_full()
            .map(|(_edge, _, _, edge_weight)| {
                if edge_weight.is_null_base() {
                    0
                } else {
                    edge_weight.copy_num
                }
            })
            .sum()
    }
    ///
    /// Set copy numbers of edges in full according to CopyNums vector
    ///
    pub fn set_copy_nums(&mut self, copy_nums: &CopyNums) {
        assert_eq!(copy_nums.len(), self.n_edges_compact());
        for edge_compact in self.compact.edge_references() {
            let copy_num = copy_nums[edge_compact.id()];

            // for all edge in compact, change copy_num of corresponding edges in full
            for &edge_in_full in edge_compact.weight().edges_in_full() {
                self.full.edge_weight_mut(edge_in_full).unwrap().copy_num = copy_num;
            }
        }
        assert!(self.is_copy_nums_valid());
    }
    ///
    /// Get current copy numbers vector
    ///
    pub fn get_copy_nums(&self) -> CopyNums {
        let mut copy_nums = CopyNums::new(self.n_edges_compact(), 0);
        for (edge, _, _, _) in self.edges_compact() {
            copy_nums[edge] = self.copy_num_of_edge_in_compact(edge);
        }
        copy_nums
    }
    ///
    /// Determine copy number of edge in compact graph
    ///
    /// = copy number of the first corresponding edge in full graph
    ///
    fn copy_num_of_edge_in_compact(&self, edge_in_compact: EdgeIndex) -> CopyNum {
        let edge_in_full = self.graph_compact()[edge_in_compact].edges_in_full()[0];
        self.graph_full()[edge_in_full].copy_num
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
/// k+1 extension
///
impl MultiDbg {
    /// Extend k to k+1.
    ///
    /// ```text
    /// k=3
    ///              TCAx2
    ///                    ┌──┐                          e1(A,2)
    ///                 ┌─►│CA│
    /// ┌──┐      ┌──┐  │  └──┘                 e0(C,3)     ┌──►v8
    /// │  │      │  ├──┘                                   │
    /// │AT├─────►│TC│               ==       v5──────►v4───┤
    /// │  │      │  ├──┐                                   │
    /// └──┘      └──┘  │  ┌──┐                             └──►v6
    ///     ATCx3       └─►│CG│
    ///                    └──┘                          e2(G,1)
    ///              TCGx1
    ///
    /// k=4
    ///         ATCAx2                               ex(A,2)
    ///               ┌───┐
    ///     ┌───┐  ┌─►│TCA│                              ┌──►v1
    ///     │   ├──┘  └───┘                              │
    ///     │ATC│                    ==           v0─────┤
    ///     │   ├──┐  ┌───┐                              │
    ///     └───┘  └─►│TCG│                              └──►v2
    ///               └───┘
    ///         ATCGx1                               ey(G,1)
    /// ```
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
    ///
    /// Upconvert path (circuit) in k MultiDbg into path in k+1 MultiDbg
    ///
    /// path is a sequence of edges in full graph.
    ///
    /// path_k corresponds to node sequence in k+1, so create edge sequence by picking a unique
    /// edge between two node.
    ///
    /// ```text
    /// k
    ///   e0    e1    e2          en-1
    ///  ----> ----> ----> ...... ---->
    ///
    /// k+1
    ///   v0 -> v1 -> v2 -> .. -> vn-1 ->
    ///      e0'   e1'   e2'  en-2'   en-1'
    /// ```
    ///
    pub fn path_kp1_from_path_k(&self, path_k: &[EdgeIndex]) -> Vec<EdgeIndex> {
        let n = path_k.len();
        let to_node_index = |edge: EdgeIndex| NodeIndex::new(edge.index());
        let mut path = Vec::new();

        for i in 0..n {
            let v_a = to_node_index(path_k[i]);
            let v_b = to_node_index(path_k[(i + 1) % n]);

            // edge from v_a to v_b
            let e = self.graph_full().find_edge(v_a, v_b).expect("invalid path");
            path.push(e);
        }

        path
    }
    /// Upconvert edge subset in k MultiDbg into edge subset in k+1 MultiDbg with corresponding
    /// emissions.
    ///
    /// edge_set_k is a set of edges in k MultiDbg.
    /// a edge in k corresponds to edges in k+1 whose target node is corresponding to the edge.
    ///
    /// Order will be preserved. If edge_set_k is sorted by some ordering, edge_set_kp1 is also
    /// sorted, although ordering in multiple corresponding edges in k+1 of an edge of k cannot be
    /// determined.
    ///
    pub fn edge_set_kp1_from_edge_set_k(&self, edge_set_k: &[EdgeIndex]) -> Vec<EdgeIndex> {
        let mut edge_set = Vec::new();
        let to_node_in_kp1 = |edge: EdgeIndex| NodeIndex::new(edge.index());

        for &edge_in_k in edge_set_k {
            let node = to_node_in_kp1(edge_in_k);
            for (edge_in_kp1, _, _) in self.childs_full(node) {
                edge_set.push(edge_in_kp1);
            }
        }

        edge_set
    }
    /// self and other is the same dbg, ignoring the node index.
    ///
    pub fn is_equivalent(&self, other: &MultiDbg) -> bool {
        self.k() == other.k()
            && self.n_nodes_full() == other.n_nodes_full()
            && self.n_edges_full() == other.n_edges_full()
            && self.n_nodes_compact() == other.n_nodes_compact()
            && self.n_edges_compact() == other.n_edges_compact()
            && self.genome_size() == other.genome_size()
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
    /// Construct node centric (full) graph G'
    ///
    /// ```text
    /// node v in G' == edge (k-mer) v in G (id will be conserved)
    /// edge from node v to node w in G' == (target of edge v is source of edge w in G)
    /// ```
    ///
    /// ```text
    /// to_node_centric_graph(G) = G'
    ///
    /// G = (V, E)
    ///
    /// G' = (V', E')
    /// V' = { f(e) | ∀e ∈ E }
    /// E' = { (f(eA), f(eB)) | ∀eA,eB ∈ E s.t. target(eA) == source(eB) }
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
    ///
    /// Purge edges whose copy numbers is 0x.
    /// This causes changes of edge index, so mapping is also returned.
    ///
    pub fn purge_edges(
        &mut self,
        edges_in_compact: &[EdgeIndex],
    ) -> (HashMap<EdgeIndex, EdgeIndex>, HashMap<EdgeIndex, EdgeIndex>) {
        // moved edges
        let mut map_full = HashMap::default();
        let mut map_compact = HashMap::default();

        // retain_edges or remove_edge
        for &edge in edges_in_compact {
            let edge_weight = self.compact.remove_edge(edge).unwrap();
            let last_edge = EdgeIndex::new(self.compact.edge_count() - 1);
            map_compact.insert(last_edge, edge);
        }

        // remove isolated nodes

        (map_full, map_compact)
    }
}

// pub fn purge_edges(graph: &mut DiGraph<N, E>)

///
/// serialize/deserialize and debug print methods
///
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
        println!("k={}", self.k());
        println!("genome_size={}", self.genome_size());
        println!("is_copy_nums_valid={}", self.is_copy_nums_valid());
        println!("degree_stats={:?}", self.degree_stats());
        println!("to_string={}", self.to_string());

        println!("Full:");
        println!("n_nodes={}", self.n_nodes_full());
        println!("n_edges={}", self.n_edges_full());
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
        println!("n_nodes={}", self.n_nodes_compact());
        println!("n_edges={}", self.n_edges_compact());
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
    /// create DBG string with `to_dbg_writer`
    ///
    pub fn to_dbg_string(&self) -> String {
        let mut writer = Vec::with_capacity(128);
        self.to_dbg_writer(&mut writer).unwrap();
        String::from_utf8(writer).unwrap()
    }
    ///
    /// create DBG file with `to_dbg_writer`
    ///
    pub fn to_dbg_file<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        let mut file = std::fs::File::create(path).unwrap();
        self.to_dbg_writer(&mut file)
    }
    /// DBG format
    ///
    /// ```text
    /// # comment
    /// # K: kmer-size
    /// K 4
    ///
    /// # N: node
    /// # id km1mer
    /// N 0  NNN
    /// N 1  ATC
    ///
    /// # E: edge
    /// # id source target seq       copy_num edges
    /// E 0  0      1      ATCGATGCT 10       8,7,2,3,1
    /// E 0  0      1      ATCGATGCT 5        8,7,2,3,1
    ///
    /// # P: path (sequence of edges)
    /// P 0,5,2,3,1,2,4,2
    ///
    /// # C: copy numbers
    /// # copy_nums_vector
    /// C 1,1,1,0,0,0,0,0,1
    /// ```
    pub fn to_dbg_writer<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        writeln!(writer, "K\t{}", self.k())?;
        for (node, weight) in self.nodes_compact() {
            writeln!(writer, "N\t{}\t{}", node.index(), self.km1mer_compact(node))?
        }
        for (edge, s, t, weight) in self.edges_compact() {
            writeln!(
                writer,
                "E\t{}\t{}\t{}\t{}\t{}\t{}",
                edge.index(),
                s.index(),
                t.index(),
                self.kmer_compact(edge),
                self.copy_num_of_edge_in_compact(edge),
                weight.edges_in_full().iter().map(|e| e.index()).format(","),
            )?
        }
        Ok(())
    }
    ///
    ///
    ///
    pub fn from_dbg_reader<R: std::io::BufRead>(reader: R) -> Self {
        let mut k = None;
        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        // sum of seq in compact graph and the number of edges in full graph
        let mut n_bases = 0;

        for line in reader.lines() {
            let text = line.unwrap();
            let first_char = text.chars().nth(0).unwrap();
            match first_char {
                'K' => {
                    let mut iter = text.split_whitespace();
                    iter.next().unwrap(); // 'K'
                    k = Some(iter.next().unwrap().parse().unwrap());
                }
                'N' => {
                    let mut iter = text.split_whitespace();
                    iter.next().unwrap(); // 'N'

                    let node: NodeIndex<DefaultIx> =
                        NodeIndex::new(iter.next().unwrap().parse().unwrap());
                    let km1mer = iter.next().unwrap().as_bytes().to_vec();

                    assert_eq!(nodes.len(), node.index(), "node is not sorted");
                    nodes.push((node, km1mer));
                }
                'E' => {
                    let k = k.unwrap();
                    let mut iter = text.split_whitespace();
                    iter.next().unwrap(); // 'E'

                    let edge: EdgeIndex<DefaultIx> =
                        EdgeIndex::new(iter.next().unwrap().parse().unwrap());
                    let s: NodeIndex<DefaultIx> =
                        NodeIndex::new(iter.next().unwrap().parse().unwrap());
                    let t: NodeIndex<DefaultIx> =
                        NodeIndex::new(iter.next().unwrap().parse().unwrap());

                    let mut kmer: Vec<u8> = iter.next().unwrap().as_bytes().to_vec();
                    let seq = kmer.split_off(k - 1);
                    n_bases += seq.len();

                    let copy_num: CopyNum = iter.next().unwrap().parse().unwrap();
                    let edges_in_full: Vec<EdgeIndex<DefaultIx>> = iter
                        .next()
                        .unwrap()
                        .split(',')
                        .map(|s| EdgeIndex::new(s.parse().unwrap()))
                        .collect();
                    assert_eq!(
                        edges_in_full.len(),
                        seq.len(),
                        "length of seq and edges_in_full is different"
                    );

                    assert_eq!(edges.len(), edge.index(), "edge is not sorted");
                    edges.push((edge, s, t, seq, copy_num, edges_in_full));
                }
                '#' => {} // pass
                _ => panic!("invalid DBG format"),
            }
        }

        // full
        let mut full = DiGraph::new();
        for (_, km1mer) in nodes.iter() {
            let is_terminal = km1mer.iter().all(|&x| x == NULL_BASE);
            full.add_node(MultiFullNode::new(is_terminal));
        }
        let mut edges_full = vec![None; n_bases];
        for (_, s, t, seq, copy_num, edges_in_full) in edges.iter() {
            // Compact:
            // s ------------------> t
            //    e0,e1,e2...
            //    ATT...
            //
            // into
            //
            // Full:
            //
            // s ----> v ----> v' ----> ... ----> v'' ----> t
            //   e0      e1       e2
            //   A       T        T
            //
            let n = seq.len();
            let mut w_prev = None;
            for i in 0..n {
                let base = seq[i];
                let edge_in_full = edges_in_full[i];

                let v = if i == 0 {
                    *s
                } else {
                    // previous w in i-1
                    w_prev.unwrap()
                };

                let w = if i == n - 1 {
                    *t
                } else {
                    // new node
                    full.add_node(MultiFullNode::new(false))
                };

                edges_full[edge_in_full.index()] =
                    Some((v, w, MultiFullEdge::new(base, *copy_num)));

                // w will be source (v) in i+1
                w_prev = Some(w);
            }
        }
        for (i, e) in edges_full.into_iter().enumerate() {
            let (source, target, weight) = e.expect("index of edge in full is wrong");
            let edge = full.add_edge(source, target, weight);
            assert_eq!(i, edge.index());
        }

        // compact
        let mut compact = DiGraph::new();
        for _ in nodes.iter() {
            compact.add_node(MultiCompactNode::new());
        }
        for (_, s, t, _, _, edges_in_full) in edges.into_iter() {
            compact.add_edge(s, t, MultiCompactEdge::new(edges_in_full));
        }

        MultiDbg {
            k: k.expect("no K section"),
            full,
            compact,
        }
    }
    ///
    /// parse DBG string with `from_dbg_reader`
    ///
    pub fn from_dbg_str(s: &str) -> Self {
        Self::from_dbg_reader(s.as_bytes())
    }
    ///
    /// parse DBG file with `from_dbg_reader`
    ///
    pub fn from_dbg_file<P: AsRef<std::path::Path>>(path: P) -> Self {
        let file = std::fs::File::open(path).unwrap();
        let reader = std::io::BufReader::new(file);
        Self::from_dbg_reader(reader)
    }
    /// GFA format
    ///
    /// ```text
    /// # comment
    ///
    /// # segment
    /// #  name  sequence   copy_number
    /// S  0     ATCGATTCG  CN:i:1
    /// S  1     ATCGATTCG  CN:i:1
    ///
    /// # link
    /// #  source  orient  target  orient  optional  node_id
    /// L  0       +       1       +       *         ID:Z:0
    /// ```
    ///
    /// * segment for each edge
    /// * link from in_edge to out_edge of each node
    ///
    pub fn to_gfa_writer<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        for (edge, s, t, weight) in self.edges_compact() {
            let seq = &self.seq_compact(edge);
            writeln!(
                writer,
                "S\t{}\t{}\tCN:i:{}\tLB:z:{}",
                edge.index(),
                sequence_to_string(&seq),
                self.copy_num_of_edge_in_compact(edge),
                sequence_to_string(&seq),
            )?
        }
        for (node, weight) in self.nodes_compact() {
            for (in_edge, _, _) in self.parents_compact(node) {
                for (out_edge, _, _) in self.childs_compact(node) {
                    writeln!(
                        writer,
                        "L\t{}\t+\t{}\t+\t*\t{}",
                        in_edge.index(),
                        out_edge.index(),
                        node.index(),
                    )?
                }
            }
        }
        Ok(())
    }
    ///
    /// create GFA string with `to_gfa_writer`
    ///
    pub fn to_gfa_string(&self) -> String {
        let mut writer = Vec::with_capacity(128);
        self.to_gfa_writer(&mut writer).unwrap();
        String::from_utf8(writer).unwrap()
    }
    ///
    /// create GFA file with `to_gfa_writer`
    ///
    pub fn to_gfa_file<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        let mut file = std::fs::File::create(path).unwrap();
        self.to_gfa_writer(&mut file)
    }
}

impl std::fmt::Display for MultiDbg {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let seqs = self
            .to_styled_seqs()
            .iter()
            .map(|seq| seq.to_string())
            .join(",");
        write!(f, "{},{}", self.k(), seqs)
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

//
// tests
//

#[cfg(test)]
mod tests {
    use super::toy;
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
    fn assert_dbg_dumpload_is_correct(dbg: MultiDbg) {
        dbg.show_graph_with_kmer();

        let s = dbg.to_dbg_string();
        println!("{}", s);

        dbg.to_dbg_file("hoge.dbg");

        let dbg_1 = MultiDbg::from_dbg_str(&s);
        dbg_1.show_graph_with_kmer();
        let s_1 = dbg.to_dbg_string();

        assert!(dbg.is_equivalent(&dbg_1));
        assert!(s == s_1);
    }
    #[test]
    fn dumpload() {
        assert_dbg_dumpload_is_correct(toy::circular());
        assert_dbg_dumpload_is_correct(toy::linear());
        assert_dbg_dumpload_is_correct(toy::intersection());
        assert_dbg_dumpload_is_correct(toy::selfloop());
        assert_dbg_dumpload_is_correct(toy::repeat());
    }
    #[test]
    fn gfa() {
        let dbg = toy::repeat();
        dbg.show_graph_with_kmer();
        let s = dbg.to_gfa_string();
        println!("{}", s);

        dbg.to_gfa_file("repeat.gfa");
    }
}
