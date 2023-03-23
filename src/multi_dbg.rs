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
    sequence_to_string, CopyNum, ReadCollection, Reads, Seq, SeqStyle, Sequence, StyledSequence,
    NULL_BASE,
};
use crate::dbg::dbg::{Dbg, DbgEdgeBase, DbgNode};
use crate::dbg::hashdbg_v2::HashDbg;
use crate::graph::compact::compact_simple_paths_for_targeted_nodes;
use crate::graph::seq_graph::{SeqEdge, SeqGraph, SeqNode};
use crate::graph::utils::{degree_stats, delete_isolated_nodes, purge_edges_with_mapping};
use crate::hmmv2::{common::PModel, hint::Hint, params::PHMMParams, table::MAX_ACTIVE_NODES};
use crate::kmer::{
    common::{kmers_to_string, KmerLike, NullableKmer},
    kmer::styled_sequence_to_kmers,
    veckmer::VecKmer,
};

use arrayvec::ArrayVec;
use fnv::FnvHashMap as HashMap;
use itertools::Itertools;
use petgraph::graph::{DefaultIx, DiGraph, EdgeIndex, NodeIndex};
use petgraph::visit::{EdgeRef, IntoNodeReferences};
use petgraph::Direction;
use petgraph_algos::iterators::{ChildEdges, EdgesIterator, NodesIterator, ParentEdges};
use rustflow::min_flow::{
    base::FlowEdgeBase, enumerate_neighboring_flows, residue::UpdateInfo, Flow,
};

pub mod draft;
pub mod output;
pub mod posterior;
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
#[derive(Clone, Debug, PartialEq)]
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
#[derive(Clone, Debug, PartialEq)]
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
#[derive(Clone, Debug, PartialEq)]
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
    ///
    /// Constructor
    ///
    pub fn new(edges_in_full: Vec<EdgeIndex>) -> Self {
        Self { edges_in_full }
    }
    ///
    /// Correspoinding edges in full graph of the edge in compact graph
    ///
    pub fn edges_in_full(&self) -> &[EdgeIndex] {
        &self.edges_in_full
    }
}

///
/// Node of MultiDbg compact graph
///
/// Empty struct
///
#[derive(Clone, Debug, PartialEq)]
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
/// # `node`
/// * compact to full: [`MultiDbg::node_in_compact_to_full`]
/// * full to compact: N/A
///     Some of nodes in full have no corresponding node in compact.
///     For terminal node, use [`MultiDbg::terminal_node_compact`]
///
/// # `terminal_node`
/// * compact: [`MultiDbg::terminal_node_compact`]
/// * full: [`MultiDbg::terminal_node_full`]
///
/// # `edge`
/// * compact to full: [`MultiCompactEdge::edges_in_full()`] stores the all corresponding edges in
/// full.
/// * full to compact
///     * single edge: [`MultiDbg::edge_in_full_to_compact`]
///     * multiple edges: [`MultiDbg::to_edge_map_into_compact`]
///
/// # `path`
/// * compact to full: [`MultiDbg::to_path_in_full`]
/// * full to compact: [`MultiDbg::to_path_in_compact`]
///
/// # copy_num
///
impl MultiDbg {
    ///
    /// Determine corresponding node in full
    ///
    pub fn node_in_compact_to_full(&self, node_in_compact: NodeIndex) -> NodeIndex {
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
    /// Get terminal node (NNNN) in compact graph
    ///
    /// The underlying genome is circular, the terminal node can be missing.
    ///
    pub fn terminal_node_compact(&self) -> Option<NodeIndex> {
        self.graph_compact().node_indices().find(|&node| {
            let v = self.node_in_compact_to_full(node);
            self.graph_full()[v].is_terminal
        })
    }
    ///
    ///
    ///
    pub fn edge_in_full_to_compact(&self, edge_in_full: EdgeIndex) -> EdgeIndex {
        let (edge, _i) = self.edge_in_full_to_compact_with_index(edge_in_full);
        edge
    }
    ///
    /// Wrapper of [`MultiCompactEdge::edges_in_full`]
    ///
    pub fn edges_in_full(&self, edge_in_compact: EdgeIndex) -> &[EdgeIndex] {
        self.graph_compact()
            .edge_weight(edge_in_compact)
            .unwrap()
            .edges_in_full()
    }
    ///
    /// Map edge in full into (edge in compact, index of edges_in_full of the edge).
    ///
    pub fn edge_in_full_to_compact_with_index(
        &self,
        edge_in_full: EdgeIndex,
    ) -> (EdgeIndex, usize) {
        self.edges_compact()
            .find_map(|(edge_in_compact, _, _, ew)| {
                ew.edges_in_full()
                    .iter()
                    .position(|&edge| edge == edge_in_full)
                    .map(|i| (edge_in_compact, i))
            })
            .expect("no corresponding edge in full")
    }
    ///
    /// Create a map from edge in full into edge in compact
    ///
    pub fn to_edge_map_into_compact(&self) -> HashMap<EdgeIndex, EdgeIndex> {
        let mut hm = HashMap::with_capacity_and_hasher(self.n_edges_full(), Default::default());
        for (edge_in_compact, _, _, ew) in self.edges_compact() {
            for edge_in_full in ew.edges_in_full() {
                hm.insert(*edge_in_full, edge_in_compact);
            }
        }
        hm
    }
}

/// Path: a sequence of edges
///
pub type Path = Vec<EdgeIndex>;

///
/// Path and styled seq conversion related
///
impl MultiDbg {
    /// Circuit (a sequence of edges in full graph) into styled sequence
    ///
    /// If the starting node of path is terminal (NNNN), it is regared as a Linear sequence.
    /// Otherwise it will be a Circular sequence.
    ///
    pub fn circuit_to_styled_seq(&self, path: &[EdgeIndex]) -> StyledSequence {
        assert!(path.len() > 0);

        let terminal_node = self.terminal_node_full();
        let mut seq = Vec::new();

        // first edge
        let (start_node, _) = self
            .graph_full()
            .edge_endpoints(path[0])
            .expect("edge in path is not in graph");
        let style = if terminal_node.is_some() && start_node == terminal_node.unwrap() {
            SeqStyle::Linear
        } else {
            SeqStyle::Circular
        };
        let mut node = start_node;

        for &e in path {
            let (s, t) = self
                .graph_full()
                .edge_endpoints(e)
                .expect("edge in path is not in graph");
            assert_eq!(
                s, node,
                "source of path[i] does not match the terminal of path[i-1]"
            );

            let w = &self.graph_full()[e];
            if w.is_null_base() {
                assert!(!style.is_circular(), "");
            } else {
                seq.push(w.base);
            }

            node = t;
        }

        assert_eq!(
            node, start_node,
            "path is not circular i.e. path[0] != path[n-1]"
        );

        StyledSequence::new(seq, style)
    }
    ///
    /// Generate Euler circuit path of full graph that traverses all edges (in full graph)
    ///
    pub fn get_euler_circuits(&self) -> Vec<Path> {
        let terminal_node = self.terminal_node_full();

        // set of paths to be returned
        let mut paths = Vec::new();
        // vector to store how many times the edge can be visited?
        let mut copy_nums_remain = self.get_copy_nums_full();

        let copy_num_of_node = |node: NodeIndex, copy_nums: &CopyNums| -> CopyNum {
            self.childs_full(node).map(|(e, _, _)| &copy_nums[e]).sum()
        };

        // pick a remaining node
        let pick_node = |copy_nums: &CopyNums| -> Option<NodeIndex> {
            // check terminal node first if remaining
            if terminal_node.is_some() && copy_num_of_node(terminal_node.unwrap(), copy_nums) > 0 {
                Some(terminal_node.unwrap())
            } else {
                self.graph_full()
                    .node_indices()
                    .find(|&v| copy_num_of_node(v, copy_nums) > 0)
            }
        };

        // pick a remaining child edge
        let pick_child =
            |node: NodeIndex, copy_nums: &CopyNums| -> Option<(EdgeIndex, NodeIndex)> {
                self.childs_full(node)
                    .sorted_by_key(|(_, _, w)| w.base)
                    .find(|(edge, _, _)| copy_nums[*edge] > 0)
                    .map(|(edge, child, _)| (edge, child))
            };

        while let Some(start_node) = pick_node(&copy_nums_remain) {
            let mut path = Vec::new();
            let mut node = start_node;

            while let Some((edge, child)) = pick_child(node, &copy_nums_remain) {
                path.push(edge);
                copy_nums_remain[edge] -= 1;
                node = child;

                if node == start_node {
                    break;
                }
            }

            if node == start_node {
                paths.push(path);
            } else {
                panic!("found path was not euler circuit");
            }
        }

        paths
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
    /// Create a set of StyledSequence that represents this MultiDbg
    ///
    /// by euler circuits that uses edges by the number of their copy numbers
    ///
    pub fn to_styled_seqs(&self) -> Vec<StyledSequence> {
        self.get_euler_circuits()
            .into_iter()
            .map(|circuit| self.circuit_to_styled_seq(&circuit))
            .collect()
    }
    /// Path validity check
    ///
    ///
    pub fn is_valid_path(&self, path_in_compact: &[EdgeIndex]) -> bool {
        let mut is_valid = true;

        for i in 0..path_in_compact.len() {
            let (_, t) = self
                .graph_compact()
                .edge_endpoints(path_in_compact[i])
                .unwrap();
            let (s, _) = self
                .graph_compact()
                .edge_endpoints(path_in_compact[(i + 1) % path_in_compact.len()])
                .unwrap();

            if t != s {
                is_valid = false;
            }
        }

        is_valid
    }
    /// Convert a path in full graph into a path in compact graph
    ///
    ///
    pub fn to_path_in_compact(&self, path_in_full: &[EdgeIndex]) -> Path {
        assert!(!path_in_full.is_empty(), "path is empty");
        let n = path_in_full.len();
        let mut path_in_compact = Vec::new();

        // detect first node
        let (_, offset) = self.edge_in_full_to_compact_with_index(path_in_full[0]);
        let mut i = 0;

        while i < n {
            let index = (n - offset + i) % n;

            // path[index]
            let edge_in_full = path_in_full[index];
            let (edge_in_compact, o) = self.edge_in_full_to_compact_with_index(edge_in_full);
            assert_eq!(o, 0, "path[index] is not the first edge in compact");
            path_in_compact.push(edge_in_compact);

            // check the following edges do match the edge_in_compact
            for (index_in_edge, &e) in self.edges_in_full(edge_in_compact).into_iter().enumerate() {
                assert_eq!(
                    path_in_full[(index + index_in_edge) % n],
                    e,
                    "ordering of edges in path is inconsistent with compact"
                );
            }

            i += self.edges_in_full(edge_in_compact).len();
        }

        assert!(self.is_valid_path(&path_in_compact));

        path_in_compact
    }
    /// Convert a path in compact graph into a path in full graph
    ///
    ///
    pub fn to_path_in_full(&self, path_in_compact: &[EdgeIndex]) -> Path {
        assert!(!path_in_compact.is_empty(), "path is empty");

        let mut path_in_full = Vec::new();
        let (mut node, _) = self
            .graph_compact()
            .edge_endpoints(path_in_compact[0])
            .unwrap();

        for &edge_in_compact in path_in_compact {
            let (s, t) = self
                .graph_compact()
                .edge_endpoints(edge_in_compact)
                .unwrap();
            assert_eq!(s, node, "invalid path in compact");

            for &edge_in_full in self.graph_compact()[edge_in_compact].edges_in_full() {
                path_in_full.push(edge_in_full);
            }

            node = t;
        }

        path_in_full
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
        let node_in_full = self.node_in_compact_to_full(node_in_compact);
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

        for &edge_in_full in self.edges_in_full(edge_in_compact).iter() {
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

        for &edge_in_full in self.edges_in_full(edge_in_compact).iter() {
            let base = self.graph_full().edge_weight(edge_in_full).unwrap().base;
            seq.push(base);
        }

        seq
    }
    ///
    /// Create kmer mapping f: Kmer -> EdgeIndex
    ///
    /// This funciton will be inefficient when k is large.
    /// If you need a path in MultiDbg of k is large, consider extending path from small k.
    ///
    pub fn to_kmer_map(&self) -> HashMap<VecKmer, EdgeIndex> {
        let mut hm = HashMap::default();
        for edge in self.graph_full().edge_indices() {
            let kmer = self.kmer_full(edge);
            hm.insert(kmer, edge);
        }
        hm
    }
    /// Convert styled seqs into paths
    ///
    pub fn paths_from_styled_seqs<T>(&self, seqs: T) -> Result<Vec<Path>, KmerNotFoundError>
    where
        T: IntoIterator,
        T::Item: AsRef<StyledSequence>,
    {
        let m = self.to_kmer_map();
        let mut paths = Vec::new();
        let mut missing_kmers = Vec::new();

        for seq in seqs {
            let mut path = Vec::new();
            for kmer in styled_sequence_to_kmers(seq.as_ref(), self.k()) {
                match m.get(&kmer) {
                    None => missing_kmers.push(kmer),
                    Some(&edge) => path.push(edge),
                }
            }
            paths.push(path);
        }

        if missing_kmers.is_empty() {
            Ok(paths)
        } else {
            Err(KmerNotFoundError(missing_kmers))
        }
    }
    /// Convert styled seqs into paths in compact graph
    ///
    pub fn compact_paths_from_styled_seqs<T>(&self, seqs: T) -> Result<Vec<Path>, KmerNotFoundError>
    where
        T: IntoIterator,
        T::Item: AsRef<StyledSequence>,
    {
        self.paths_from_styled_seqs(seqs).map(|paths| {
            paths
                .into_iter()
                .map(|path| self.to_path_in_compact(&path))
                .collect()
        })
    }
}

#[derive(Clone, Debug)]
pub struct KmerNotFoundError(Vec<VecKmer>);

impl std::fmt::Display for KmerNotFoundError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "KmerNotFoundError({})", kmers_to_string(&self.0))
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
        assert!(self.is_copy_nums_valid(), "invalid new copy_nums");
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
    /// Get neighboring copy numbers
    ///
    ///
    pub fn to_neighbor_copy_nums_and_infos(
        &self,
        max_cycle_size: usize,
        max_flip: usize,
    ) -> Vec<(CopyNums, UpdateInfo)> {
        let network = self.graph_compact().map(
            |_, _| (),
            |edge, _| {
                let copy_num = self.copy_num_of_edge_in_compact(edge);
                FlowEdgeBase::new(copy_num.saturating_sub(1), copy_num.saturating_add(1), 0.0)
            },
        );
        let copy_nums = self.get_copy_nums();
        enumerate_neighboring_flows(&network, &copy_nums, Some(max_cycle_size), Some(max_flip))
    }
    ///
    ///
    pub fn to_neighbor_copy_nums(&self, max_cycle_size: usize, max_flip: usize) -> Vec<CopyNums> {
        self.to_neighbor_copy_nums_and_infos(max_cycle_size, max_flip)
            .into_iter()
            .map(|(copy_nums, _)| copy_nums)
            .collect()
    }
    ///
    /// Paths in compact Vec<Path> into CopyNums
    ///
    pub fn copy_nums_from_compact_path<P: AsRef<[EdgeIndex]>, PS: AsRef<[P]>>(
        &self,
        paths: PS,
    ) -> CopyNums {
        let mut copy_nums = CopyNums::new(self.n_edges_compact(), 0);

        for path in paths.as_ref() {
            assert!(self.is_valid_path(path.as_ref()));

            for &edge in path.as_ref() {
                copy_nums[edge] += 1;
            }
        }

        copy_nums
    }
}

///
/// k+1 extension
///
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
            true,
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
    /// if e[0] == nnnA  and e[n-1] == Gnnn, a terminal node v == nnnn should be added to path in
    /// order for e to be valid path in k+1.
    /// ```text
    /// k=4             Linear                          k=4            Circular
    ///
    ///         en-2    en-1    e0      e1                      e0      e1      e2
    ///    ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐           ┌───┐   ┌───┐   ┌───┐   ┌───┐
    /// ───▶AGn├───▶Gnn├───▶nnn├───▶nnA├───▶nAA├───▶       │AGG├───▶GGA├───▶GAA├───▶AAG│
    ///    └───┘   └───┘   └───┘   └───┘   └───┘           └─▲─┘   └───┘   └───┘   └─┬─┘
    ///                                                      └───────────────────────┘
    /// k=5                                             k=5             e3
    ///         vn-2    vn-1    v0      v1                      v0      v1      v2
    ///       ┌────┐  ┌────┐  ┌────┐  ┌────┐                  ┌────┐  ┌────┐  ┌────┐
    ///    ───▶AGnn├──▶Gnnn│  │nnnA├──▶nnAA├──▶               │AGGA├──▶GGAA├──▶GAAG│
    ///       └────┘  └─┬──┘  └──▲─┘  └────┘                  └─▲──┘  └────┘  └──┬─┘
    ///                 │        │                              │     ┌────┐     │
    ///                 │        │                              └─────┤AAGG◀─────┘
    ///                 │ ┌────┐ │                                    └────┘
    ///                 └─▶nnnn├─┘                                      v3
    ///                   └────┘
    ///                     v
    /// ```
    ///
    pub fn path_kp1_from_path_k(&self, path_k_in_full: &[EdgeIndex]) -> Path {
        let n = path_k_in_full.len();
        let to_node_index = |edge: EdgeIndex| NodeIndex::new(edge.index());
        let mut path = Vec::new();

        // check path[0] is starting node and path[n-1] is ending node
        // starting node
        //   = nnnA
        //   = parent is terminal (nnn)
        //
        // ending node
        //   = Gnnn
        //   = child is terminal (nnn)
        let terminal = self.terminal_node_full();
        let first = to_node_index(path_k_in_full[0]);
        let last = to_node_index(path_k_in_full[n - 1]);
        let start = terminal.and_then(|t| self.graph_full().find_edge(t, first));
        let end = terminal.and_then(|t| self.graph_full().find_edge(last, t));

        if start.is_some() && end.is_some() {
            // lienar
            path.push(start.unwrap());
            for i in 0..(n - 1) {
                let v_a = to_node_index(path_k_in_full[i]);
                let v_b = to_node_index(path_k_in_full[i + 1]);

                // edge from v_a to v_b
                let e = self.graph_full().find_edge(v_a, v_b).expect("invalid path");
                path.push(e);
            }
            path.push(end.unwrap());
        } else {
            // circular
            for i in 0..n {
                let v_a = to_node_index(path_k_in_full[i]);
                let v_b = to_node_index(path_k_in_full[(i + 1) % n]);

                // edge from v_a to v_b
                let e = self.graph_full().find_edge(v_a, v_b).expect("invalid path");
                path.push(e);
            }
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
    /// ```text
    ///             k=4              k=5
    ///
    ///             │
    ///           ┌─▼─┐               │
    ///           │ATC│      ┌─▶XATCC │
    ///           └─┬─┘      │      ┌─▼──┐
    ///             │  ATCC──┘      │ATCC│
    ///           ┌─▼─┐             └─┬──┘
    ///           │TCC│      ┌─▶ATCCG │
    ///           └─┬─┘      │      ┌─▼──┐
    ///             │  TCCG──┘      │TCCG│
    ///           ┌─▼─┐             └─┬──┘
    ///           │CCG│   base        │
    ///           └─┬─┘   corresponds ▼
    ///             │
    ///             ▼
    /// ```
    ///
    pub fn edge_set_kp1_from_edge_set_k(&self, edge_set_k: &[EdgeIndex]) -> Vec<EdgeIndex> {
        let mut edge_set = Vec::new();
        let to_node_in_kp1 = |edge: EdgeIndex| NodeIndex::new(edge.index());

        for &edge_in_k in edge_set_k {
            let node = to_node_in_kp1(edge_in_k);
            for (edge_in_kp1, _, _) in self.parents_full(node) {
                println!("node={:?} edge={:?}", node, edge_in_kp1);
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
    /// exactly same graph including node/edge index
    ///
    pub fn is_equal(&self, other: &MultiDbg) -> bool {
        // full
        let sfn: Vec<_> = self.nodes_full().collect();
        let ofn: Vec<_> = other.nodes_full().collect();
        let sfe: Vec<_> = self.edges_full().collect();
        let ofe: Vec<_> = other.edges_full().collect();
        // compact
        let scn: Vec<_> = self.nodes_compact().collect();
        let ocn: Vec<_> = other.nodes_compact().collect();
        let sce: Vec<_> = self.edges_compact().collect();
        let oce: Vec<_> = other.edges_compact().collect();
        sfn == ofn && sfe == ofe && scn == ocn && sce == oce
    }
}

///
/// Profile HMM related
///
impl MultiDbg {
    ///
    ///
    ///
    fn to_seq_graph(&self) -> DiGraph<SNode, SEdge> {
        self.to_node_centric_graph(
            |e, ew| SNode {
                copy_num: ew.copy_num,
                base: ew.base,
            },
            || {
                let copy_num_terminal = self.copy_num_of_node(
                    self.terminal_node_full()
                        .expect("there is no terminal node"),
                );
                SNode {
                    copy_num: copy_num_terminal,
                    base: NULL_BASE,
                }
            },
            |_, _, _| SEdge {},
            |_| SEdge {},
            false,
        )
    }
    /// Convert MultiDbg into node-centric Profile HMM `PModel` via SeqGraph `DiGraph<SNode, SEdge>`
    ///
    ///
    pub fn to_phmm(&self, param: PHMMParams) -> PModel {
        self.to_seq_graph().to_phmm(param)
    }
    ///
    ///
    pub fn to_uniform_phmm(&self, param: PHMMParams) -> PModel {
        self.to_seq_graph().to_uniform_phmm(param)
    }
}

///
/// Minimum instance of SeqNode from MultiDbg
///
struct SNode {
    copy_num: CopyNum,
    base: u8,
}

///
/// Minimum instance of SeqEdge from MultiDbg
///
struct SEdge {}

impl SeqNode for SNode {
    fn copy_num(&self) -> CopyNum {
        self.copy_num
    }
    fn base(&self) -> u8 {
        self.base
    }
}

impl SeqEdge for SEdge {
    fn copy_num(&self) -> Option<CopyNum> {
        None
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
    /// * `to_node(e: edge_index, edge_weight)`
    ///     For all edge (k-mer) e in G, `to_node(e)` is called and the k-mer is added to G' as
    ///     node.
    /// * `to_terminal_node()`
    ///     A special k-mer `NNNN` will be newly added to G' while no corresponding edge in G.
    /// * `to_edge(e+: edge_index_source, e-: edge_index_target, v: center_node_index)`
    ///     For all pair of in/out edges (e+, e-) of all node v in G, edge from e+ to e- will be
    ///     added to G'.
    /// * `to_terminal_edge(e: edge)`
    ///     For all pair of in/out edges of terminal node NNN in G, edge will be added to G'
    ///     corresponding to each edge in G.
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
        add_terminal: bool,
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
                if add_terminal {
                    let terminal_node = graph.add_node(to_terminal_node());
                    for (e, _, _) in self.parents_full(node) {
                        let v = to_node_index(e);
                        graph.add_edge(v, terminal_node, to_terminal_edge(e));
                    }
                    for (e, _, _) in self.childs_full(node) {
                        let v = to_node_index(e);
                        graph.add_edge(terminal_node, v, to_terminal_edge(e));
                    }
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
    /// `(map_full, map_compact)`
    ///
    /// `map_{full,compact}` is a mapping from edge in original graph into edge in modified graph or None if
    /// deleted.
    ///
    pub fn purge_edges(
        &mut self,
        edges_in_compact: &[EdgeIndex],
    ) -> (
        HashMap<EdgeIndex, Option<EdgeIndex>>,
        HashMap<EdgeIndex, Option<EdgeIndex>>,
    ) {
        // List up edges to be removed in full
        let mut edges_in_full = Vec::new();
        for &edge in edges_in_compact {
            for &edge_in_full in self.graph_compact()[edge].edges_in_full() {
                edges_in_full.push(edge_in_full);
            }
        }

        // remove edges from full/compact graph
        let (map_compact, _) = purge_edges_with_mapping(&mut self.compact, edges_in_compact);
        let (map_full, _) = purge_edges_with_mapping(&mut self.full, &edges_in_full);

        // remove isolated nodes
        delete_isolated_nodes(&mut self.compact);
        delete_isolated_nodes(&mut self.full);

        // update edges_in_full in compact::MultiCompactEdge
        for weight in self.compact.edge_weights_mut() {
            for e in weight.edges_in_full.iter_mut() {
                *e = match map_full.get(&e) {
                    None => *e,
                    Some(&v) => v.expect("remove edge twice"),
                };
            }
        }

        (map_full, map_compact)
    }
    ///
    ///
    ///
    pub fn purge_and_extend(
        &self,
        edges_in_compact_to_purge: &[EdgeIndex],
        paths: Option<Vec<Path>>,
        hints: Option<Vec<Hint>>,
    ) -> (Self, Option<Vec<Path>>, Option<Vec<Hint>>) {
        // Delete edges from the graph k-DBG and get mappings between edges
        //
        let mut dbg_k = self.clone();
        let (map_full, _) = dbg_k.purge_edges(&edges_in_compact_to_purge);

        // (2)
        // Extend to k+1-DBG
        //
        let dbg_kp1 = dbg_k.to_kp1_dbg();

        // (3)
        // convert paths and hints
        //
        // (a) paths
        let paths = match paths {
            Some(paths_k_compact) => paths_k_compact
                .into_iter()
                .map(|path_k_compact| {
                    let path_k_full_purged: Option<Path> = self
                        .to_path_in_full(&path_k_compact)
                        .into_iter()
                        .map(|e| map_full.get(&e).copied().unwrap_or(Some(e)))
                        .collect();

                    match path_k_full_purged {
                        None => None,
                        Some(path_k_compact_purged) => {
                            let path_kp1_full =
                                dbg_kp1.path_kp1_from_path_k(&path_k_compact_purged);
                            let path_kp1_compact = dbg_kp1.to_path_in_compact(&path_kp1_full);
                            Some(path_kp1_compact)
                        }
                    }
                })
                .collect(),
            None => None,
        };

        // (b) hints
        //
        let hints = match hints {
            Some(hints) => Some(
                hints
                    .into_iter()
                    .map(|hint_k| {
                        // hint is Vec<ArrayVec<NodeIndex>>
                        //
                        // NodeIndex corresponds to EdgeIndex in k
                        // 1. index-change caused by purge
                        //
                        // 2. index-change caused by extension
                        //
                        let vecs = hint_k.to_inner();
                        Hint::new(
                            vecs.into_iter()
                                .map(|vec| {
                                    // 1
                                    // nodes in hmm -> edges in purged k dbg
                                    let edges_in_k: Vec<EdgeIndex> = vec
                                        .into_iter()
                                        .filter_map(|node_in_hmm| {
                                            // map purged
                                            let edge_in_full = EdgeIndex::new(node_in_hmm.index());
                                            let edge_in_purged = map_full
                                                .get(&edge_in_full)
                                                .copied()
                                                .unwrap_or(Some(edge_in_full));
                                            edge_in_purged
                                        })
                                        .collect();

                                    // 2
                                    // edges in k -> edges in k+1 -> nodes in hmm
                                    // into arrayvec by taking top n elements
                                    dbg_kp1
                                        .edge_set_kp1_from_edge_set_k(&edges_in_k)
                                        .into_iter()
                                        .map(|e| NodeIndex::new(e.index()))
                                        .take(MAX_ACTIVE_NODES)
                                        .collect::<ArrayVec<NodeIndex, MAX_ACTIVE_NODES>>()
                                })
                                .collect(),
                        )
                    })
                    .collect(),
            ),
            None => None,
        };

        (dbg_kp1, paths, hints)
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::toy;
    use super::*;
    use crate::common::{ei, ni};
    use crate::dbg::mocks as dbgmocks;
    use crate::kmer::veckmer::kmer;

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
    #[test]
    fn purge_edges_for_toy() {
        let mut dbg = toy::intersection();
        dbg.show_graph_with_kmer();
        assert_eq!(dbg.kmer_compact(ei(0)), kmer(b"ATCAnnn"));
        assert_eq!(dbg.kmer_compact(ei(1)), kmer(b"ATCCnnn"));
        assert_eq!(dbg.kmer_compact(ei(2)), kmer(b"nnnTATC"));
        assert_eq!(dbg.kmer_compact(ei(3)), kmer(b"nnnGATC"));

        let (mf, mc) = dbg.purge_edges(&[ei(1), ei(2)]);
        dbg.show_graph_with_kmer();
        println!("{:?} {:?}", mf, mc);
        assert_eq!(dbg.n_edges_full(), 8);
        assert_eq!(dbg.n_edges_compact(), 2);
        assert_eq!(dbg.kmer_compact(ei(0)), kmer(b"ATCAnnn"));
        assert_eq!(dbg.kmer_compact(ei(1)), kmer(b"nnnGATC"));
    }
    #[test]
    fn neighbors_for_toy() {
        let mut dbg = toy::intersection();
        dbg.show_graph_with_kmer();
        let neighbors = dbg.to_neighbor_copy_nums_and_infos(10, 0);
        for (copy_nums, update_info) in neighbors {
            println!("{} {:?}", copy_nums, update_info);
        }

        {
            let neighbors = dbg.to_neighbor_copy_nums(10, 0);
            assert_eq!(neighbors.len(), 8);
        }

        {
            let neighbors = dbg.to_neighbor_copy_nums(10, 2);
            assert_eq!(neighbors.len(), 12);
        }
    }
    #[test]
    fn phmm_for_toy() {
        let dbg = toy::intersection();
        dbg.show_graph_with_kmer();

        let phmm = dbg.to_phmm(PHMMParams::uniform(0.01));
        println!("{}", phmm);
    }
    #[test]
    fn phmm_from_dbg() {
        let param = PHMMParams::uniform(0.01);
        let reads = &[b"CTAGCTT"];
        let dbg = dbgmocks::mock_intersection();
        let phmm_dbg = dbg.to_phmm(param);
        println!("{}", phmm_dbg);
        let p_dbg = phmm_dbg.to_full_prob(reads);
        println!("{}", p_dbg);

        let mdbg: MultiDbg = dbg.into();
        let phmm_mdbg = mdbg.to_phmm(param);
        println!("{}", phmm_mdbg);
        let p_mdbg = phmm_mdbg.to_full_prob(reads);
        println!("{}", p_mdbg);

        assert_eq!(p_dbg, p_mdbg);
    }
    #[test]
    fn bridge_full_and_compact() {
        {
            let dbg = toy::circular();
            dbg.show_graph_with_kmer();
            assert_eq!(dbg.terminal_node_full(), None);
            assert_eq!(dbg.terminal_node_compact(), None);
            assert_eq!(dbg.edge_in_full_to_compact(ei(2)), ei(0));
            assert_eq!(dbg.edge_in_full_to_compact_with_index(ei(2)), (ei(0), 1));

            // path
            let p = dbg.to_path_in_full(&[ei(0)]);
            assert_eq!(p, vec![ei(1), ei(2), ei(3), ei(0)]);

            let q = dbg.to_path_in_compact(&p);
            println!("{:?}", q);
            assert_eq!(q, vec![ei(0)]);

            let q = dbg.to_path_in_compact(&vec![ei(2), ei(3), ei(0), ei(1)]);
            println!("{:?}", q);
            assert_eq!(q, vec![ei(0)]);
        }

        {
            let dbg = toy::repeat();
            dbg.show_graph_with_kmer();
            // nnn
            assert_eq!(dbg.terminal_node_full(), Some(ni(0)));
            assert_eq!(dbg.terminal_node_compact(), Some(ni(0)));
            // CAG
            assert_eq!(dbg.node_in_compact_to_full(ni(1)), ni(6));
            // multiple edge
            let m = dbg.to_edge_map_into_compact();
            println!("{:?}", m);
            for (&edge_in_full, &edge_in_compact) in m.iter() {
                assert!(dbg.graph_compact()[edge_in_compact]
                    .edges_in_full()
                    .contains(&edge_in_full));
            }
            // single edge
            assert_eq!(dbg.edge_in_full_to_compact(ei(14)), ei(0));
            assert_eq!(dbg.edge_in_full_to_compact_with_index(ei(3)), (ei(2), 3));

            // path
            println!("#1");
            let p = dbg.to_path_in_full(&[ei(2), ei(0)]);
            println!("{:?}", p);
            assert_eq!(
                p,
                vec![
                    ei(0),
                    ei(1),
                    ei(2),
                    ei(3),
                    ei(4),
                    ei(5),
                    ei(9),
                    ei(10),
                    ei(11),
                    ei(12),
                    ei(13),
                    ei(14)
                ]
            );
            let q = dbg.to_path_in_compact(&p);
            println!("{:?}", q);
            assert_eq!(q, vec![ei(2), ei(0)]);

            // path starts from non-terminal edge
            println!("#2");
            let p = vec![
                ei(2),
                ei(3),
                ei(4),
                ei(5),
                ei(9), // p[4]
                ei(10),
                ei(11),
                ei(12),
                ei(13),
                ei(14),
                ei(0), // p[10]
                ei(1),
            ];
            let q = dbg.to_path_in_compact(&p);
            println!("{:?}", q);
            assert_eq!(q, vec![ei(2), ei(0)]);
        }
    }
    #[test]
    fn styled_seq() {
        {
            let dbg = toy::circular();
            dbg.show_graph_with_kmer();

            let c = dbg.get_euler_circuits();
            println!("{:?}", c);
            assert_eq!(c, vec![vec![ei(0), ei(1), ei(2), ei(3)]]);

            let s = dbg.to_styled_seqs();
            println!("{:?}", s);
            assert_eq!(s, vec![StyledSequence::circular(b"GATC".to_vec())]);

            let p = dbg.paths_from_styled_seqs(&s);
            assert_eq!(p.unwrap(), vec![vec![ei(3), ei(0), ei(1), ei(2)]]);
        }

        {
            let dbg = toy::repeat();
            dbg.show_graph_with_kmer();

            let c = dbg.get_euler_circuits();
            println!("{:?}", c);
            assert_eq!(
                c,
                vec![vec![
                    ei(0),
                    ei(1),
                    ei(2),
                    ei(3),
                    ei(4),
                    ei(5),
                    ei(6),
                    ei(7),
                    ei(8),
                    ei(6),
                    ei(7),
                    ei(8),
                    ei(6),
                    ei(7),
                    ei(8),
                    ei(9),
                    ei(10),
                    ei(11),
                    ei(12),
                    ei(13),
                    ei(14)
                ]]
            );

            let s = dbg.to_styled_seqs();
            println!("{:?}", s);
            assert_eq!(
                s,
                vec![StyledSequence::linear(b"TCCCAGCAGCAGCAGGAA".to_vec())]
            );

            let p = dbg.paths_from_styled_seqs(&s);
            println!("{:?}", p);
            assert_eq!(p.unwrap(), c);

            let p = dbg.paths_from_styled_seqs(vec![StyledSequence::linear(b"TCCCAGGAA".to_vec())]);
            println!("{:?}", p);
            assert_eq!(
                p.unwrap(),
                vec![vec![
                    ei(0),
                    ei(1),
                    ei(2),
                    ei(3),
                    ei(4),
                    ei(5),
                    ei(9),
                    ei(10),
                    ei(11),
                    ei(12),
                    ei(13),
                    ei(14)
                ]]
            );
        }
    }
    #[test]
    fn purge_and_extend_toy() {
        // remove an edge
        {
            let mut dbg_k = toy::repeat();
            println!("### k");
            dbg_k.set_copy_nums(&vec![1, 0, 1].into());
            dbg_k.show_graph_with_kmer();
            assert_eq!(dbg_k.n_nodes_full(), 14);
            assert_eq!(dbg_k.n_edges_full(), 15);
            assert_eq!(dbg_k.n_nodes_compact(), 2);
            assert_eq!(dbg_k.n_edges_compact(), 3);
            assert_eq!(dbg_k.genome_size(), 9);

            // purge repetitive edges
            let hints = vec![Hint::from(vec![
                vec![
                    ni(2), // nTCC -> nnTCC ni(3)
                    ni(3), // TCCC -> nTCCC ni(4)
                ],
                vec![
                    ni(3), // TCCC -> nTCCC ni(4)
                    ni(4), // CCCA -> TCCCA ni(5)
                ],
                vec![
                    ni(6), // CAGC -> none
                    ni(9), // CAGG -> CCAGG ni(7)
                ],
            ])];
            let (dbg_kp1, paths, hints) =
                dbg_k.purge_and_extend(&[ei(1)], Some(vec![vec![ei(2), ei(0)]]), Some(hints));
            println!("### k+1");
            dbg_kp1.show_graph_with_kmer();
            assert_eq!(dbg_kp1.n_nodes_full(), 13);
            assert_eq!(dbg_kp1.n_edges_full(), 13);
            assert_eq!(dbg_kp1.n_nodes_compact(), 1);
            assert_eq!(dbg_kp1.n_edges_compact(), 1);
            assert_eq!(dbg_k.genome_size(), 9);
            assert!(dbg_kp1.is_copy_nums_valid());
            println!("paths={:?}", paths);
            println!("hints={:?}", hints);
            assert_eq!(paths, Some(vec![vec![ei(0)]]));
            assert_eq!(
                hints,
                Some(vec![Hint::from(vec![
                    vec![ni(3), ni(4)],
                    vec![ni(4), ni(5)],
                    vec![ni(7)],
                ])])
            );
        }

        // no purge of edges
        {
            let dbg_k = toy::repeat();
            println!("### k");
            dbg_k.show_graph_with_kmer();
            let paths = vec![vec![ei(2), ei(1), ei(1), ei(1), ei(0)]];
            let hints = vec![Hint::from(vec![
                vec![
                    ni(2), // nTCC -> nnTCC ni(3)
                    ni(3), // TCCC -> nTCCC ni(4)
                ],
                vec![
                    ni(3), // TCCC -> nTCCC ni(4)
                    ni(4), // CCCA -> TCCCA ni(5)
                ],
                vec![
                    ni(6), // CAGC -> CCAGC ni(10) or GCAGC ni(8)
                    ni(9), // CAGG -> CCAGG ni(9)  or GCAGG ni(7)
                ],
            ])];

            let (dbg_kp1, paths, hints) = dbg_k.purge_and_extend(&[], Some(paths), Some(hints));
            println!("### k+1");
            dbg_kp1.show_graph_with_kmer();
            println!("paths={:?}", paths);
            println!("hints={:?}", hints);
            assert!(dbg_kp1.is_copy_nums_valid());
            assert!(dbg_kp1.is_equal(&toy::repeat_kp1()));
            assert_eq!(
                paths,
                Some(vec![vec![
                    ei(6),
                    ei(1),
                    ei(3),
                    ei(5),
                    ei(3),
                    ei(5),
                    ei(3),
                    ei(0),
                    ei(4)
                ]])
            );
            assert_eq!(
                hints,
                Some(vec![Hint::from(vec![
                    vec![ni(3), ni(4)],
                    vec![ni(4), ni(5)],
                    vec![ni(10), ni(8), ni(9), ni(7)],
                ])])
            );
        }
    }
}
