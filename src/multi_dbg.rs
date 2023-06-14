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
use crate::graph::euler::euler_circuit_count;
use crate::graph::k_shortest::{k_shortest_cycle, k_shortest_simple_path};
use crate::graph::seq_graph::{SeqEdge, SeqGraph, SeqNode};
use crate::graph::utils::{
    bridge_edges, degree_stats, delete_isolated_nodes, purge_edges_with_mapping, EdgeMap,
};
use crate::hmmv2::{
    common::PModel,
    hint::{Mapping, Mappings},
    params::PHMMParams,
    table::MAX_ACTIVE_NODES,
};
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
use petgraph_algos::common::is_edge_simple;
use petgraph_algos::iterators::{ChildEdges, EdgesIterator, NodesIterator, ParentEdges};
use rustflow::min_flow::{
    base::{FlowEdgeBase, FlowGraph},
    enumerate_neighboring_flows, find_neighboring_flow_by_edge_change,
    residue::{
        flow_to_residue_convex, is_meaningful_move_on_residue_graph, residue_graph_cycle_to_flow,
        ResidueDirection, UpdateInfo,
    },
    Flow,
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
#[derive(Clone, Debug, PartialEq)]
pub struct MultiCompactNode {
    ///
    /// if this node corresponds to k-1mer NNN, then true.
    ///
    is_terminal: bool,
}

impl MultiCompactNode {
    pub fn new(is_terminal: bool) -> Self {
        Self { is_terminal }
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
    ///
    /// # of edges with base in full graph
    ///
    pub fn n_emittable_edges(&self) -> usize {
        self.edges_full()
            .filter(|(_, _, _, ew)| !ew.is_null_base())
            .count()
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
    /// Get terminal node (= (k-1)-mer NNN) in full graph
    ///
    /// The underlying genome is circular, the terminal node can be missing.
    ///
    pub fn terminal_node_full(&self) -> Option<NodeIndex> {
        self.nodes_full()
            .find(|(_, node_weight)| node_weight.is_terminal)
            .map(|(node, _)| node)
    }
    ///
    /// Get terminal node (= (k-1)-mer NNN) in compact graph
    ///
    /// The underlying genome is circular, the terminal node can be missing.
    ///
    pub fn terminal_node_compact(&self) -> Option<NodeIndex> {
        self.nodes_compact()
            .find(|(_, node_weight)| node_weight.is_terminal)
            .map(|(node, _)| node)
    }
    ///
    /// Check if an edge in compact graph is either
    /// * starting edge (source node is terminal NNN)
    /// * ending edge (target node is terminal NNN)
    ///
    pub fn is_start_or_end_edge_compact(&self, edge_in_compact: EdgeIndex) -> bool {
        match self.terminal_node_compact() {
            Some(terminal) => {
                let (s, t) = self
                    .graph_compact()
                    .edge_endpoints(edge_in_compact)
                    .unwrap();
                s == terminal || t == terminal
            }
            None => false,
        }
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
    /// the number of bases the edge in compact represents?
    /// (= the number of collapsed edges in full)
    ///
    pub fn n_bases(&self, edge_in_compact: EdgeIndex) -> usize {
        self.edges_in_full(edge_in_compact).len()
    }
    ///
    /// Number of linear haplotypes
    /// = sum of in/out-degree of the terminal node
    ///
    pub fn n_haplotypes(&self) -> usize {
        match self.terminal_node_full() {
            Some(terminal_node) => self.copy_num_of_node(terminal_node),
            None => 0,
        }
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
    ///
    /// Dump fasta of Eulerian traverse
    ///
    pub fn to_fasta<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bio::io::fasta::Writer::new(file);

        for (i, s) in self.to_styled_seqs().iter().enumerate() {
            writer.write(&format!("p{}", i), Some(&s.style().to_string()), s.seq())?;
        }

        Ok(())
    }
    /// Path (in compact graph) validity check
    ///
    pub fn is_valid_compact_path(&self, path_in_compact: &[EdgeIndex]) -> bool {
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
    /// Path (in full graph) validity check
    ///
    pub fn is_valid_full_path(&self, path_in_full: &[EdgeIndex]) -> bool {
        let mut is_valid = true;

        for i in 0..path_in_full.len() {
            let (_, t) = self.graph_full().edge_endpoints(path_in_full[i]).unwrap();
            let (s, _) = self
                .graph_full()
                .edge_endpoints(path_in_full[(i + 1) % path_in_full.len()])
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

        assert!(self.is_valid_compact_path(&path_in_compact));

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
    ///
    /// Paths in full Vec<Path> into CopyNums
    ///
    pub fn copy_nums_from_full_path<P: AsRef<[EdgeIndex]>, PS: AsRef<[P]>>(
        &self,
        paths: PS,
    ) -> CopyNums {
        let mut copy_nums = CopyNums::new(self.n_edges_compact(), 0);

        for path in paths.as_ref() {
            let compact_path = self.to_path_in_compact(path.as_ref());

            for edge in compact_path {
                copy_nums[edge] += 1;
            }
        }

        copy_nums
    }
    ///
    /// Calculate log #EC = (logarithm of number of Euler circuits)
    ///
    /// see [Kingsford2010](https://doi.org/10.1186/1471-2105-11-21)
    ///
    pub fn n_euler_circuits(&self) -> f64 {
        let graph = self
            .graph_compact()
            .map(|_, _| (), |e, _| self.copy_num_of_edge_in_compact(e));

        euler_circuit_count(&graph)
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
    /// Get maximum copy number of edges (k-mers)
    ///
    pub fn max_copy_num(&self) -> CopyNum {
        self.edges_full()
            .map(|(_, _, _, edge_weight)| edge_weight.copy_num)
            .max()
            .unwrap()
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
    pub fn copy_num_of_edge_in_compact(&self, edge_in_compact: EdgeIndex) -> CopyNum {
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
        // TODO use ArrayVec?
        let copy_num_ins: Vec<CopyNum> = self
            .parents_full(node)
            .map(|(_, _, w)| w.copy_num)
            .collect();
        let copy_num_outs: Vec<CopyNum> =
            self.childs_full(node).map(|(_, _, w)| w.copy_num).collect();
        let index_in = self
            .parents_full(node)
            .position(|(edge, _, _)| edge == edge_in)
            .expect("edge_in is not a parent of node");
        let index_out = self
            .childs_full(node)
            .position(|(edge, _, _)| edge == edge_out)
            .expect("edge_out is not a child of node");
        Self::guess_copy_num(&copy_num_ins, &copy_num_outs)[index_in][index_out]
    }
    /// Assign copy number to the pair of edges incoming/outgoing from the same node.
    /// (used in `guess_copy_num_of_kp1_edge`.)
    ///
    /// Give the same amount of copy number to all >0x outgoing edge in order of copy numbers of
    /// outgoing edge.
    ///
    pub fn guess_copy_num(
        copy_num_ins: &[CopyNum],
        copy_num_outs: &[CopyNum],
    ) -> Vec<Vec<CopyNum>> {
        let n_in = copy_num_ins.len();
        let n_out = copy_num_outs.len();
        assert_eq!(
            copy_num_ins.iter().sum::<CopyNum>(),
            copy_num_outs.iter().sum::<CopyNum>()
        );

        // copy_nums[index_in][index_out] is the copy number of k+1 edge (edge_in -> edge_out)
        let mut copy_nums = vec![vec![0; n_out]; n_in];
        let mut remain_ins = copy_num_ins.to_owned();
        let mut remain_outs = copy_num_outs.to_owned();

        while remain_ins.iter().any(|&x| x > 0) && remain_outs.iter().any(|&x| x > 0) {
            for i_in in 0..n_in {
                for i_out in 0..n_out {
                    if remain_ins[i_in] > 0 && remain_outs[i_out] > 0 {
                        copy_nums[i_in][i_out] += 1;
                        remain_ins[i_in] -= 1;
                        remain_outs[i_out] -= 1;
                    }
                }
            }
        }

        copy_nums
    }
    /// Get neighboring copy numbers
    ///
    ///
    pub fn to_neighbor_copy_nums_and_infos(
        &self,
        config: NeighborConfig,
    ) -> Vec<(CopyNums, UpdateInfo)> {
        // (1) enumerate all short cycles
        let mut short_cycles = self.to_short_neighbors(config.max_cycle_size, config.max_flip);
        eprintln!("short_cycles: {}", short_cycles.len());

        // (2) add long cycles for 0x -> 1x
        if config.use_long_cycles {
            let mut long_cycles = self.to_long_neighbors();
            eprintln!("long_cycles: {}", long_cycles.len());
            short_cycles.append(&mut long_cycles);
        }

        // (3) add long cycles for reducing high-copy edges
        if config.use_reducers {
            let mut reducers = self.to_reducer_neighbors();
            eprintln!("reducers: {}", reducers.len());
            short_cycles.append(&mut reducers);
        }

        short_cycles
    }
    ///
    /// Create flow network `FlowGraph` for enumerating neighboring copy numbers (flow)
    ///
    pub fn to_flow_network(&self) -> FlowGraph<usize> {
        // create flow network for rustflow
        let network = self.graph_compact().map(
            |_, _| (),
            |edge, _| {
                let copy_num = self.copy_num_of_edge_in_compact(edge);
                FlowEdgeBase::new(copy_num.saturating_sub(1), copy_num.saturating_add(1), 0.0)
            },
        );
        network
    }
    /// Rescue neighbors
    /// 0x
    pub fn to_rescue_neighbors(
        &self,
        k: usize,
        not_make_new_zero_edge: bool,
    ) -> Vec<(CopyNums, UpdateInfo)> {
        let mut ret = vec![];
        for (edge, _, _, _) in self.edges_compact() {
            if self.copy_num_of_edge_in_compact(edge) == 0 {
                for n in self.to_rescue_neighbors_for_edge(edge, k, not_make_new_zero_edge) {
                    ret.push(n);
                }
            }
        }
        ret
    }
    /// Rescue neighbors
    ///
    /// 0x
    pub fn to_rescue_neighbors_for_edge(
        &self,
        edge: EdgeIndex,
        k: usize,
        not_make_new_zero_edge: bool,
    ) -> Vec<(CopyNums, UpdateInfo)> {
        let copy_nums = self.get_copy_nums();
        let network = self.graph_compact().map(
            |_, _| (),
            |edge, _| {
                let copy_num = self.copy_num_of_edge_in_compact(edge);
                if copy_num == 0 {
                    FlowEdgeBase::new(0, copy_num.saturating_add(1), 0.0)
                } else if not_make_new_zero_edge {
                    FlowEdgeBase::new(1, copy_num.saturating_add(1), 0.0)
                } else {
                    FlowEdgeBase::new(copy_num.saturating_sub(1), copy_num.saturating_add(1), 0.0)
                }
            },
        );
        // println!("{:?}", petgraph::dot::Dot::with_config(&network, &[]));
        let rg = flow_to_residue_convex(&network, &copy_nums);
        let edge_in_rg = rg
            .edge_indices()
            .find(|&e| rg[e].target == edge && rg[e].direction == ResidueDirection::Up)
            .unwrap();
        let (v, w) = rg.edge_endpoints(edge_in_rg).unwrap();
        let cycles = k_shortest_simple_path(&rg, w, v, k, |e| {
            if e == edge_in_rg {
                usize::MAX
            } else {
                self.n_bases(rg[e].target)
            }
        });

        cycles
            .into_iter()
            .map(|mut cycle| {
                cycle.insert(0, edge_in_rg);
                cycle
            })
            .filter(|cycle| is_edge_simple(&rg, &cycle))
            .map(|cycle| residue_graph_cycle_to_flow(&copy_nums, &rg, &cycle))
            .collect()
    }
    ///
    ///
    ///
    pub fn to_short_neighbors(
        &self,
        max_cycle_size: usize,
        max_flip: usize,
    ) -> Vec<(CopyNums, UpdateInfo)> {
        let network = self.to_flow_network();
        let copy_nums = self.get_copy_nums();

        enumerate_neighboring_flows(&network, &copy_nums, Some(max_cycle_size), Some(max_flip))
    }
    ///
    ///
    ///
    pub fn to_long_neighbors(&self) -> Vec<(CopyNums, UpdateInfo)> {
        let network = self.to_flow_network();
        let copy_nums = self.get_copy_nums();

        // let max_copy_num = self.max_copy_num();
        self.graph_compact()
            .edge_indices()
            .filter(|&e| self.copy_num_of_edge_in_compact(e) == 0)
            .filter_map(|e| {
                find_neighboring_flow_by_edge_change(
                    &network,
                    &copy_nums,
                    e,
                    ResidueDirection::Up,
                    // weight function
                    |e| self.n_bases(e) / (self.copy_num_of_edge_in_compact(e) + 1),
                    // |e| max_copy_num - self.copy_num_of_edge_in_compact(e),
                )
            })
            // .filter(|(_copy_nums, info)| info.len() > config.max_cycle_size)
            .filter(|(_copy_nums, info)| !self.is_passing_terminal(&info))
            .collect()
    }
    pub fn to_reducer_neighbors(&self) -> Vec<(CopyNums, UpdateInfo)> {
        let network = self.graph_compact().map(
            |_, _| (),
            |edge, _| {
                let copy_num = self.copy_num_of_edge_in_compact(edge);
                if copy_num > 2 {
                    // reduceable
                    FlowEdgeBase::new(copy_num.saturating_sub(1), copy_num, 0.0)
                } else {
                    FlowEdgeBase::new(copy_num, copy_num, 0.0)
                }
            },
        );
        let copy_nums = self.get_copy_nums();
        enumerate_neighboring_flows(&network, &copy_nums, Some(100), Some(0))
    }
    pub fn is_passing_terminal(&self, info: &UpdateInfo) -> bool {
        info.iter()
            .any(|&(edge, _)| self.is_start_or_end_edge_compact(edge))
    }
    pub fn has_zero_to_one_change(&self, info: &UpdateInfo) -> bool {
        info.iter().any(|&(edge, dir)| {
            self.copy_num_of_edge_in_compact(edge) == 0 && dir == ResidueDirection::Up
        })
    }
}

///
/// Parameters for neighbor search of copy numbers
///
#[derive(Clone, Debug, Copy)]
pub struct NeighborConfig {
    /// Max size of cycle (in compact graph) in BFS short-cycle search
    ///
    pub max_cycle_size: usize,
    /// Max number of flips (+/- or -/+ changes) in cycles on compact residue graph
    ///
    pub max_flip: usize,
    /// Augment short cycles with long cycles causing 0x -> 1x change
    ///
    pub use_long_cycles: bool,
    /// Ignore cyclic paths passing through the terminal node NNN
    /// because this changes the number of haplotypes.
    ///
    pub ignore_cycles_passing_terminal: bool,
    ///
    /// Long cycle with Down direction only to reduce high copy number edges
    ///
    pub use_reducers: bool,
}

///
/// k+1 extension and mappings
///
/// # Extend into k+1
/// * to_kp1_dbg
/// * path_kp1_from_path_k
/// * edge_set_kp1_from_edge_set_k
/// * hint_kp1_from_hint_k
///
/// # Extend k unless no ambiguity
/// * to_k_max_dbg
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
    ///
    /// Upconvert Mapping (edge subset) in k MultiDbg into edge subset in k+1 MultiDbg with corresponding
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
    /// node in k-HMM == edge in k-HMM == node in k+1 DBG
    ///
    pub fn hint_kp1_from_hint_k(&self, hint_k: &Mapping) -> Mapping {
        hint_k.map_nodes(|node_in_k_hmm| {
            let node_in_kp1 = node_in_k_hmm;
            let mut nodes_in_kp1_hmm = Vec::new();
            for (edge_in_kp1, _, _) in self.parents_full(node_in_kp1) {
                nodes_in_kp1_hmm.push(NodeIndex::new(edge_in_kp1.index()));
            }
            nodes_in_kp1_hmm
        })
    }
    ///
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
            |_node, node_weight| MultiCompactNode::new(node_weight.is_terminal),
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
}

///
/// Edge purging and mappings
///
impl MultiDbg {
    ///
    /// Purge edges whose copy numbers is 0x.
    /// This causes changes of edge index, so mapping is also returned.
    ///
    /// Note that purging edges in compact can make simple-path collapsable nodes in the
    /// resulting compact graph, but this function does not perform re-simple-path-collapsing.
    ///
    /// Return type is `PurgeEdgeMap`
    /// Wrapper of `HashMap<EdgeIndex, Option<EdgeIndex>>`
    /// that maps edge in graph before purging into edge in graph after purging.
    ///
    pub fn purge_edges(&mut self, edges_in_compact: &[EdgeIndex]) -> PurgeEdgeMap {
        //
        // [1] update compact graph
        //
        // remove edges from compact graph
        let mut map_compact = purge_edges_with_mapping(&mut self.compact, edges_in_compact);
        // removing edges may cause isolated edges (i.e. bridging edges)
        // so remove them next.
        let bridge_edges = bridge_edges(&self.compact);
        // println!("bridge={:?}", bridge_edges);
        // let map_compact_bridge = purge_edges_with_mapping(&mut self.compact, &bridge_edges);
        // map_compact.compose(&map_compact_bridge);

        //
        // [2] update full graph
        //
        // a. List up edges to be removed in full graph, corresponding edges_in_compact and bridge_edges
        // in compact graph.
        let mut edges_in_full = Vec::new();
        for &edge in edges_in_compact.iter().chain(bridge_edges.iter()) {
            for &edge_in_full in self.graph_compact()[edge].edges_in_full() {
                edges_in_full.push(edge_in_full);
            }
        }
        // b. remove edges from full graph
        let map_full = purge_edges_with_mapping(&mut self.full, &edges_in_full);

        // remove isolated nodes
        delete_isolated_nodes(&mut self.compact);
        delete_isolated_nodes(&mut self.full);

        //
        // [3] update links between compact and full
        //
        // update edges_in_full in compact::MultiCompactEdge
        for weight in self.compact.edge_weights_mut() {
            // renew old indexes in edges_in_full weight of each edge
            for e in weight.edges_in_full.iter_mut() {
                *e = map_full
                    .from_original(*e)
                    .unwrap_or_else(|| panic!("e{} is deleted", e.index()));
            }
        }

        PurgeEdgeMap {
            map_full,
            map_compact,
        }
    }
    ///
    /// Purge edges and extend DBG into k=k_max.
    ///
    /// With index conversion of paths (that represents genome) and hints (for efficient phmm of
    /// reads) caused by purging and extending.
    ///
    /// * paths: set of path (in full graph)
    /// * hints: set of hint (node list of phmm, that corresponds to edges in full)
    ///
    /// If stop_when_ambiguous=true, stop extending when ambiguous (i.e. some nodes have
    /// in_deg/out_deg > 1, so copy numbers of k+1 DBG is not unique) even if k has not reached k_max.
    ///
    pub fn purge_and_extend(
        &self,
        edges_in_compact_to_purge: &[EdgeIndex],
        k_max: usize,
        stop_when_ambiguous: bool,
        paths: Option<Vec<Path>>,
        mappings: &Mappings,
    ) -> (Self, Option<Vec<Path>>, Mappings) {
        let mut dbg_k = self.clone();

        //
        // Delete edges from the graph k-DBG and get mappings between edges
        //
        let m = dbg_k.purge_edges(&edges_in_compact_to_purge);
        let mut paths: Option<Vec<Path>> =
            paths.and_then(|paths| paths.into_iter().map(|path| m.update_path(&path)).collect());
        let mut mappings = Mappings::new(
            mappings
                .into_iter()
                .map(|hint| m.update_mapping(&hint))
                .collect(),
        );

        // (2) Extend DBG
        //
        while dbg_k.k() < k_max {
            // Extend
            //
            assert!(dbg_k.is_copy_nums_valid());
            let dbg_kp1 = dbg_k.to_kp1_dbg();
            assert!(dbg_kp1.is_copy_nums_valid());
            // Convert paths and hints
            //
            // (a) paths
            paths = paths.map(|paths| {
                paths
                    .into_iter()
                    .map(|path| dbg_kp1.path_kp1_from_path_k(&path))
                    .collect()
            });
            // (b) mappings
            mappings = Mappings::new(
                mappings
                    .into_iter()
                    .map(|hint_k| dbg_kp1.hint_kp1_from_hint_k(hint_k))
                    .collect(),
            );

            //
            let is_ambiguous = dbg_k.n_ambiguous_node() > 0;
            dbg_k = dbg_kp1;
            if stop_when_ambiguous && is_ambiguous {
                break;
            }
        }

        (dbg_k, paths, mappings)
    }
}

/// Mapping of edges in full/compact before and after edge-purging.
///
/// Return type of [`MultiDbg::purge_edges`]
///
/// the map (`map_compact` and `map_full`) can be generated by `purge_edges_with_mapping`.
///
#[derive(Clone, Debug)]
pub struct PurgeEdgeMap {
    map_compact: EdgeMap,
    map_full: EdgeMap,
}

impl PurgeEdgeMap {
    ///
    /// Map compact edge in G into compact edge in G' (after edge purging)
    ///
    pub fn compact(&self, edge_in_compact: EdgeIndex) -> Option<EdgeIndex> {
        self.map_compact.from_original(edge_in_compact)
    }
    ///
    /// Map full edge in G into full edge in G' (after edge purging)
    ///
    pub fn full(&self, edge_in_full: EdgeIndex) -> Option<EdgeIndex> {
        self.map_full.from_original(edge_in_full)
    }
    ///
    /// Convert Path = Vec<EdgeIndex> using edge-map before-and-after edge purging.
    ///
    pub fn update_path(&self, path_in_full: &[EdgeIndex]) -> Option<Vec<EdgeIndex>> {
        path_in_full.iter().map(|&e| self.full(e)).collect()
    }
    ///
    /// Convert Mapping
    ///
    pub fn update_mapping(&self, mapping: &Mapping) -> Mapping {
        mapping.map_nodes(|node| {
            let edge = EdgeIndex::new(node.index());
            match self.full(edge) {
                Some(edge_after) => vec![NodeIndex::new(edge_after.index())],
                None => vec![],
            }
        })
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
    use crate::e2e::{generate_dataset, ReadType};
    use crate::genome;
    use crate::kmer::veckmer::kmer;
    use crate::prob::p;

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
    fn copy_num_guess() {
        assert_eq!(
            MultiDbg::guess_copy_num(&[2], &[1, 0, 1]),
            vec![vec![1, 0, 1]]
        );

        assert_eq!(
            MultiDbg::guess_copy_num(&[2, 6, 0, 4], &[12, 0]),
            vec![vec![2, 0], vec![6, 0], vec![0, 0], vec![4, 0]]
        );

        assert_eq!(
            MultiDbg::guess_copy_num(&[2, 6, 0, 4], &[9, 0, 1, 2]),
            vec![
                vec![1, 0, 1, 0],
                vec![5, 0, 0, 1],
                vec![0, 0, 0, 0],
                vec![3, 0, 0, 1],
            ]
        );
    }
    #[test]
    fn purge_edges_for_toy() {
        let mut dbg = toy::intersection();
        dbg.show_graph_with_kmer();
        assert_eq!(dbg.kmer_compact(ei(0)), kmer(b"ATCAnnn"));
        assert_eq!(dbg.kmer_compact(ei(1)), kmer(b"ATCCnnn"));
        assert_eq!(dbg.kmer_compact(ei(2)), kmer(b"nnnTATC"));
        assert_eq!(dbg.kmer_compact(ei(3)), kmer(b"nnnGATC"));

        let m = dbg.purge_edges(&[ei(1), ei(2)]);
        dbg.show_graph_with_kmer();
        println!("{:?}", m);
        assert_eq!(m.compact(ei(1)), None);
        assert_eq!(m.compact(ei(2)), None);
        assert_eq!(dbg.n_edges_full(), 8);
        assert_eq!(dbg.n_edges_compact(), 2);
        assert_eq!(dbg.kmer_compact(ei(0)), kmer(b"ATCAnnn"));
        assert_eq!(dbg.kmer_compact(ei(1)), kmer(b"nnnGATC"));
    }
    #[test]
    fn neighbors_for_toy() {
        let mut dbg = toy::intersection();
        dbg.show_graph_with_kmer();
        let neighbors = dbg.to_neighbor_copy_nums_and_infos(NeighborConfig {
            max_cycle_size: 10,
            max_flip: 0,
            use_long_cycles: false,
            ignore_cycles_passing_terminal: false,
            use_reducers: false,
        });
        for (copy_nums, update_info) in neighbors {
            println!("{} {:?}", copy_nums, update_info);
        }

        {
            let neighbors = dbg.to_neighbor_copy_nums_and_infos(NeighborConfig {
                max_cycle_size: 10,
                max_flip: 0,
                use_long_cycles: false,
                ignore_cycles_passing_terminal: false,
                use_reducers: false,
            });
            assert_eq!(neighbors.len(), 8);
        }

        {
            let neighbors = dbg.to_neighbor_copy_nums_and_infos(NeighborConfig {
                max_cycle_size: 10,
                max_flip: 2,
                use_long_cycles: false,
                ignore_cycles_passing_terminal: false,
                use_reducers: false,
            });
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
    fn purge_and_extend_toy_with_purge() {
        // extend only once
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
            let paths = vec![dbg_k.to_path_in_full(&[ei(2), ei(0)])];
            let mappings = Mappings::new(vec![Mapping::from_nodes_and_probs(&vec![
                vec![
                    (ni(2), p(0.7)), // nTCC -> nnTCC ni(3)
                    (ni(3), p(0.3)), // TCCC -> nTCCC ni(4)
                ],
                vec![
                    (ni(3), p(0.8)), // TCCC -> nTCCC ni(4)
                    (ni(4), p(0.2)), // CCCA -> TCCCA ni(5)
                ],
                vec![
                    (ni(6), p(0.6)), // CAGC -> none
                    (ni(9), p(0.4)), // CAGG -> CCAGG ni(7)
                ],
            ])]);
            println!("{:?} {:?}", paths, mappings);

            {
                let (dbg_kp1, paths, mappings) =
                    dbg_k.purge_and_extend(&[ei(1)], 5, false, Some(paths.clone()), &mappings);
                println!("### k+1");
                dbg_kp1.show_graph_with_kmer();
                assert_eq!(dbg_kp1.k(), 5);
                assert_eq!(dbg_kp1.n_nodes_full(), 13);
                assert_eq!(dbg_kp1.n_edges_full(), 13);
                assert_eq!(dbg_kp1.n_nodes_compact(), 1);
                assert_eq!(dbg_kp1.n_edges_compact(), 1);
                assert_eq!(dbg_kp1.genome_size(), 9);
                assert!(dbg_kp1.is_copy_nums_valid());
                println!("paths={:?}", paths);
                println!("mappings={:?}", mappings);
                assert_eq!(paths, Some(vec![dbg_kp1.to_path_in_full(&[ei(0)])]));
                assert_eq!(
                    mappings,
                    Mappings::new(vec![Mapping::from_nodes_and_probs(&vec![
                        vec![(ni(3), p(0.7)), (ni(4), p(0.3))],
                        vec![(ni(4), p(0.8)), (ni(5), p(0.2))],
                        vec![(ni(7), p(0.4))],
                    ])]),
                );
            }

            {
                println!("### k=10");
                let (dbg_k10, paths_k10, mappings_k10) = dbg_k.purge_and_extend(
                    &[ei(1)],
                    10,
                    true, // there is no ambiguity from k=4 to k=10, so k=k_max can be reached
                    Some(paths.clone()),
                    &mappings,
                );
                dbg_k10.show_graph_with_kmer();
                println!("paths={:?}", paths_k10);
                println!("mappings={:?}", mappings_k10);

                // check dbg
                assert_eq!(dbg_k10.k(), 10);
                assert_eq!(dbg_k10.genome_size(), 9);
                assert_eq!(dbg_k10.n_nodes_full(), 18);
                assert_eq!(dbg_k10.n_edges_full(), 18);
                assert_eq!(dbg_k10.n_nodes_compact(), 1);
                assert_eq!(dbg_k10.n_edges_compact(), 1);
                assert_eq!(
                    dbg_k10.kmer_compact(ei(0)),
                    VecKmer::from_bases(b"nnnnnnnnnTCCCAGGAAnnnnnnnnn")
                );
                assert_eq!(
                    dbg_k10.km1mer_compact(ni(0)),
                    VecKmer::from_bases(b"nnnnnnnnn")
                );
                assert!(dbg_k10.is_copy_nums_valid());
                // check path and hint
                assert_eq!(paths_k10, Some(vec![dbg_k10.to_path_in_full(&[ei(0)])]));
                assert_eq!(
                    mappings_k10,
                    Mappings::new(vec![Mapping::from_nodes_and_probs(&vec![
                        // 10-mers ending with..
                        vec![(ni(13), p(0.7)), (ni(0), p(0.3))], // nnTCC and nTCCC
                        vec![(ni(0), p(0.8)), (ni(1), p(0.2))],  // nTCCC and TCCCA
                        vec![(ni(3), p(0.4))],                   // CCAGG
                    ])])
                );
            }
        }
    }

    #[test]
    fn purge_and_extend_toy_without_purge() {
        // no purge of edges
        {
            let dbg_k = toy::repeat();
            println!("### k");
            dbg_k.show_graph_with_kmer();
            let paths = vec![dbg_k.to_path_in_full(&[ei(2), ei(1), ei(1), ei(1), ei(0)])];
            let genome = vec![StyledSequence::linear(b"TCCCAGCAGCAGCAGGAA".to_vec())];
            assert_eq!(paths, dbg_k.paths_from_styled_seqs(&genome).unwrap());
            let mappings = Mappings::new(vec![Mapping::from_nodes_and_probs(&vec![
                vec![
                    (ni(2), p(0.7)), // nTCC -> nnTCC ni(3)
                    (ni(3), p(0.3)), // TCCC -> nTCCC ni(4)
                ],
                vec![
                    (ni(3), p(0.8)), // TCCC -> nTCCC ni(4)
                    (ni(4), p(0.2)), // CCCA -> TCCCA ni(5)
                ],
                vec![
                    (ni(6), p(0.6)), // CAGC -> CCAGC ni(10) or GCAGC ni(8)
                    (ni(9), p(0.4)), // CAGG -> CCAGG ni(9)  or GCAGG ni(7)
                ],
            ])]);

            let dbg_k5 = toy::repeat_kp1();
            let paths_k5 = vec![dbg_k5.to_path_in_full(&[
                ei(6),
                ei(1),
                ei(3),
                ei(5),
                ei(3),
                ei(5),
                ei(3),
                ei(0),
                ei(4),
            ])];
            let mappings_k5 = Mappings::new(vec![Mapping::from_nodes_and_probs(&vec![
                vec![(ni(3), p(0.7)), (ni(4), p(0.3))],
                vec![(ni(4), p(0.8)), (ni(5), p(0.2))],
                vec![
                    (ni(10), p(0.6) / 2),
                    (ni(8), p(0.6) / 2),
                    (ni(9), p(0.4) / 2),
                    (ni(7), p(0.4) / 2),
                ],
            ])]);

            {
                println!("### k+1");
                let (dbg_kp1, paths_kp1, mappings_kp1) =
                    dbg_k.purge_and_extend(&[], 5, false, Some(paths.clone()), &mappings);
                dbg_kp1.show_graph_with_kmer();
                println!("paths={:?}", paths_kp1);
                println!("mappings={:?}", mappings_kp1);

                assert_eq!(dbg_kp1.k(), 5);
                assert!(dbg_kp1.is_copy_nums_valid());
                assert!(dbg_kp1.is_equal(&toy::repeat_kp1()));
                assert_eq!(paths_kp1.unwrap(), paths_k5);
                assert_eq!(mappings_kp1, mappings_k5);

                // there is ambiguity when kp1=5 so output is same even if k_max=10.
                let (dbg_kp1, paths_kp1, mappings_kp1) =
                    dbg_k.purge_and_extend(&[], 10, true, Some(paths.clone()), &mappings);
                assert_eq!(dbg_kp1.k(), 5);
                assert!(dbg_kp1.is_copy_nums_valid());
                assert!(dbg_kp1.is_equal(&toy::repeat_kp1()));
                assert_eq!(paths_kp1.unwrap(), paths_k5);
                assert_eq!(mappings_kp1, mappings_k5);
            }

            {
                println!("### k=10");
                let (dbg_k10, paths_k10, mappings_k10) =
                    dbg_k.purge_and_extend(&[], 10, false, Some(paths.clone()), &mappings);
                dbg_k10.show_graph_with_kmer();
                println!("paths={:?}", paths_k10);
                println!("mappings={:?}", mappings_k10);

                assert_eq!(dbg_k10.k(), 10);
                assert_eq!(dbg_k10.genome_size(), 18);
                assert!(dbg_k10.is_copy_nums_valid());
                assert_eq!(
                    paths_k10.unwrap(),
                    dbg_k10.paths_from_styled_seqs(&genome).unwrap()
                );
                for m in mappings_k10.clone().into_iter() {
                    println!("m={}", m);
                }
                assert_eq!(
                    mappings_k10[0].nodes(0),
                    vec![ni(31), ni(0)], // nTCC and TCCC
                );
                assert_eq!(
                    mappings_k10[0].nodes(1),
                    vec![ni(0), ni(1)], // TCCC and CCCA
                );
                assert_eq!(
                    mappings_k10[0].nodes(2),
                    vec![ni(3), ni(4), ni(11), ni(9), ni(12), ni(10)] // CAGG and CAGC
                );
            }
        }
    }
    #[test]
    fn n_euler_circuits_test_toy() {
        {
            let dbg = toy::repeat();
            dbg.show_graph_with_kmer();
            let n = dbg.n_euler_circuits();
            println!("n={}", n.exp());
            assert!((n.exp() - 1.0).abs() < 0.0001);
        }

        {
            let mut dbg = toy::one_in_n_repeat();
            dbg.show_graph_with_kmer();
            let n = dbg.n_euler_circuits();
            println!("n={}", n.exp());
            assert!((n.exp() - 5.0).abs() < 0.0001);

            dbg.set_copy_nums(&CopyNums::new(dbg.n_edges_compact(), 0));
            dbg.show_graph_with_kmer();
            let n = dbg.n_euler_circuits();
            println!("n={}", n.exp());
            assert!((n.exp() - 0.0).abs() < 0.0001);
        }

        {
            let dbg = toy::two_components();
            dbg.show_graph_with_kmer();
            let n = dbg.n_euler_circuits();
            println!("n={}", n.exp());
            // FIXME separate for connected components
        }
    }
    #[test]
    #[ignore = "consume memory"]
    fn n_euler_circuits_test_large() {
        let genome = genome::tandem_repeat_polyploid_with_unique_homo_ends(
            100, 10, 0, 0.0, 0, 300, 2, 0.01, 0,
        );
        let param = PHMMParams::uniform(0.01);
        let dataset = generate_dataset(genome, 0, 20, 500, ReadType::FragmentWithRevComp, param);
        let dbg = MultiDbg::create_draft_from_dataset(20, &dataset);

        let n = dbg.n_euler_circuits();
        println!("n={}", n);
        assert!((n - 61.6).abs() < 0.1);
    }
}
