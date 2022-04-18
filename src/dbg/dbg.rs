//!
//! De bruijn graph definitions
//!
//!
use super::edge_centric::{
    SimpleEDbg, SimpleEDbgEdge, SimpleEDbgEdgeWithAttr, SimpleEDbgNode, SimpleEDbgWithAttr,
};
use super::impls::{SimpleDbg, SimpleDbgEdge, SimpleDbgNode};
use super::intersections::Intersection;
use crate::common::{CopyNum, Reads, SeqStyle, Sequence, StyledSequence};
use crate::dbg::hashdbg_v2::HashDbg;
use crate::graph::iterators::{ChildEdges, EdgesIterator, NodesIterator, ParentEdges};
use crate::kmer::kmer::styled_sequence_to_kmers;
use crate::kmer::{KmerLike, NullableKmer};
use crate::vector::{DenseStorage, EdgeVec, NodeVec};
use fnv::{FnvHashMap as HashMap, FnvHashSet as HashSet};
use itertools::{iproduct, izip, Itertools};
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph::Direction;

pub type NodeCopyNums = NodeVec<DenseStorage<CopyNum>>;
pub type EdgeCopyNums = EdgeVec<DenseStorage<CopyNum>>;

///
/// (Node-centric) De bruijn graph struct
/// k
///
#[derive(Clone)]
pub struct Dbg<N: DbgNode, E: DbgEdge> {
    ///
    /// k-mer size
    ///
    k: usize,
    ///
    /// backend DiGraph with node-centric representation
    ///
    pub graph: DiGraph<N, E>,
}

///
/// Trait for nodes in Dbg
///
pub trait DbgNode: Clone {
    type Kmer: KmerLike + NullableKmer;
    fn new(kmer: Self::Kmer, copy_num: CopyNum) -> Self;
    ///
    /// Kmer of this node of the Dbg
    fn kmer(&self) -> &Self::Kmer;
    ///
    /// Copy number count of this node in Dbg
    fn copy_num(&self) -> CopyNum;
    ///
    /// Modify copy number count of this node
    fn set_copy_num(&mut self, copy_num: CopyNum);
    ///
    /// Single base assigned to this node in Dbg
    /// Last base of kmer will be used as an emission
    fn emission(&self) -> u8 {
        self.kmer().last()
    }
    ///
    /// this node is emittable or not?
    /// i.e. emission is not b'N'.
    ///
    fn is_emittable(&self) -> bool {
        self.emission() != b'N'
    }
    ///
    /// check if k-mer of this node is head (NNNNA)
    ///
    fn is_head(&self) -> bool {
        self.kmer().is_head()
    }
    ///
    /// check if k-mer of this node is tail (ANNNN)
    ///
    fn is_tail(&self) -> bool {
        self.kmer().is_tail()
    }
    ///
    /// calculate the genome size of this node
    ///
    fn genome_size(&self) -> CopyNum {
        if self.is_emittable() {
            self.copy_num()
        } else {
            0
        }
    }
}

///
/// Trait for edges in Dbg
///
pub trait DbgEdge: Clone {
    fn new(copy_num: Option<CopyNum>) -> Self;
    ///
    /// Copy number count of this edge in Dbg
    /// None if number is not assigned. (i.e. random transition)
    fn copy_num(&self) -> Option<CopyNum>;
    ///
    /// Modify copy number count of this edge
    fn set_copy_num(&mut self, copy_num: Option<CopyNum>);
}

///
/// Basic graph operations for Dbg
///
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    /// k-mer size of the de Bruijn Graph
    pub fn k(&self) -> usize {
        self.k
    }
    /// create iterator of all nodes
    /// Item of the iterator is `(NodeIndex, &N)`.
    pub fn nodes(&self) -> NodesIterator<N> {
        NodesIterator::new(&self.graph)
    }
    /// create iterator of all nodes
    /// Item of the iterator is
    /// `(EdgeIndex, NodeIndex of source, NodeIndex of target, &E)`.
    pub fn edges(&self) -> EdgesIterator<E> {
        EdgesIterator::new(&self.graph)
    }
    /// create iterator of all child edges of the node
    ///
    /// Item of the iterator is `(EdgeIndex, NodeIndex, Edge weight)`
    ///
    /// * `EdgeIndex` is index of edge (node, child)
    /// * `NodeIndex` is index of child
    /// * `EdgeWeight` is of the edge transition
    pub fn childs(&self, node: NodeIndex) -> ChildEdges<E> {
        ChildEdges::new(&self.graph, node)
    }
    /// create iterator of all parent edges of the node
    ///
    /// Item of the iterator is `(EdgeIndex, NodeIndex, Edge weight)`
    ///
    /// * `EdgeIndex` is index of edge (parent, node)
    /// * `NodeIndex` is index of parent
    /// * `EdgeWeight` is of the edge transition
    pub fn parents(&self, node: NodeIndex) -> ParentEdges<E> {
        ParentEdges::new(&self.graph, node)
    }
    ///
    /// NodeIndex(s) of heads/tails
    ///
    pub fn tips(&self) -> Intersection<N::Kmer> {
        let mut in_nodes = Vec::new();
        let mut out_nodes = Vec::new();

        for (node, weight) in self.nodes() {
            // ANNN
            if weight.kmer().is_tail() {
                in_nodes.push(node);
            }
            // NNNA
            if weight.kmer().is_head() {
                out_nodes.push(node);
            }
        }

        Intersection::new(N::Kmer::null_kmer(self.k()), in_nodes, out_nodes)
    }
    /// Return the number of nodes in the graph
    pub fn n_nodes(&self) -> usize {
        self.graph.node_count()
    }
    /// Return the number of edges in the graph
    pub fn n_edges(&self) -> usize {
        self.graph.edge_count()
    }
    /// Get a reference of node weight
    pub fn node(&self, node: NodeIndex) -> &N {
        self.graph.node_weight(node).unwrap()
    }
    /// Get a reference of edge weight
    pub fn edge(&self, edge: EdgeIndex) -> &E {
        self.graph.edge_weight(edge).unwrap()
    }
    /// convert node to kmer reference
    pub fn kmer(&self, node: NodeIndex) -> &N::Kmer {
        self.graph.node_weight(node).unwrap().kmer()
    }
    /// convert node to emission
    pub fn emission(&self, node: NodeIndex) -> u8 {
        self.graph.node_weight(node).unwrap().emission()
    }
    /// check if two nodes `a, b: NodeIndex` is connected or not
    pub fn contains_edge(&self, a: NodeIndex, b: NodeIndex) -> bool {
        self.graph.contains_edge(a, b)
    }
    /// Lookup an edge from `a: NodeIndex` to `b: NodeIndex`.
    pub fn find_edge(&self, a: NodeIndex, b: NodeIndex) -> Option<EdgeIndex> {
        self.graph.find_edge(a, b)
    }
    /// The number of edges from `a: NodeIndex` to `b: NodeIndex`.
    pub fn count_edge(&self, a: NodeIndex, b: NodeIndex) -> usize {
        self.graph.edges_connecting(a, b).count()
    }
    /// Check if the edge is a warping edge or not.
    /// warping edge is edges like `ANNNN -> NNNNC`
    pub fn is_warp_edge(&self, edge: EdgeIndex) -> bool {
        let (v, w) = self.graph.edge_endpoints(edge).unwrap();
        self.kmer(v).is_tail() && self.kmer(w).is_head()
    }
    /// Lookup a node from k-mer. It takes O(V).
    pub fn find_node_from_kmer(&self, kmer: &N::Kmer) -> Option<NodeIndex> {
        self.graph.node_indices().find(|&v| self.kmer(v) == kmer)
    }
    /// Lookup an edge from k+1-mer. It takes O(E).
    pub fn find_edge_from_kp1mer(&self, kp1mer: &N::Kmer) -> Option<EdgeIndex> {
        let prefix = kp1mer.prefix();
        let suffix = kp1mer.suffix();
        self.edges()
            .find(|(_, v, w, _)| self.kmer(*v) == &prefix && self.kmer(*w) == &suffix)
            .map(|(e, _, _, _)| e)
    }
}

///
/// Basic properties
///
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    /// Get the total genome size of this de bruijn graph
    /// It can be calculated by the sum of genome sizes
    /// of all emittable nodes.
    ///
    pub fn genome_size(&self) -> CopyNum {
        self.nodes().map(|(_, weight)| weight.genome_size()).sum()
    }
    /// CopyNums of nodes are consistent, that is
    /// 'sum of copynums of childs' == 'sum of copynums of siblings'
    /// for each kmers
    pub fn has_consistent_node_copy_nums(&self) -> bool {
        self.to_edbg().has_consistent_copy_nums()
    }
    /// CopyNums of edges are consistent?
    /// i.e. the sum of copy numbers on in-edges and out-edges
    /// are equal to the copy numbers of the node.
    ///
    /// Head/tail nodes will break this consistency
    /// because warping edges (from tail to head) do not have copy number.
    pub fn has_consistent_edge_copy_nums(&self) -> bool {
        self.nodes().all(|(node, weight)| {
            let in_copy_nums: Option<CopyNum> = self
                .graph
                .edges_directed(node, Direction::Incoming)
                .map(|e| e.weight().copy_num())
                .sum();
            let out_copy_nums: Option<CopyNum> = self
                .graph
                .edges_directed(node, Direction::Outgoing)
                .map(|e| e.weight().copy_num())
                .sum();

            // if the node is head or tail, allow the inconsistency.
            if weight.is_head() || weight.is_tail() {
                true
            } else {
                in_copy_nums == out_copy_nums
            }
        })
    }
    /// Check if all the edges has a copy number
    /// except the warp edges (ANNN -> NNNC)
    ///
    pub fn is_edge_copy_nums_assigned(&self) -> bool {
        self.edges()
            .all(|(e, _, _, ew)| self.is_warp_edge(e) || ew.copy_num().is_some())
    }
    /// Check if there is no node duplicates.
    pub fn has_no_duplicated_node(&self) -> bool {
        let mut s = HashSet::default();
        for (_, weight) in self.nodes() {
            let kmer = weight.kmer().clone();
            if s.contains(&kmer) {
                return false;
            } else {
                s.insert(kmer);
            }
        }
        true
    }
    ///
    /// Check if there is no parallel edges.
    ///
    pub fn has_no_parallel_edge(&self) -> bool {
        self.edges().all(|(_, v, w, _)| self.count_edge(v, w) == 1)
    }
    ///
    /// Check if backend-graph is valid.
    ///
    /// Complexity: O(|V|^2)
    ///
    pub fn is_graph_valid(&self) -> bool {
        let m = self.to_kmer_map();
        for (v, weight) in self.nodes() {
            // if graph has child kmer, the graph should have an edge
            // from the node to the child.
            for child in weight.kmer().childs() {
                if let Some(&w) = m.get(&child) {
                    if !self.contains_edge(v, w) {
                        println!("is_graph_valid: two kmers {}/{} (node index {}/{}) do not have an edge!", weight.kmer(), child, v.index(), w.index());
                        return false;
                    }
                }
            }

            // if graph has parent kmer, the graph should have an edge
            // from the parent to the node.
            for parent in weight.kmer().parents() {
                if let Some(&w) = m.get(&parent) {
                    if !self.contains_edge(w, v) {
                        println!("is_graph_valid: two kmers {}/{} (node index {}/{}) do not have an edge!", parent, weight.kmer(), w.index(), v.index());
                        return false;
                    }
                }
            }
        }
        true
    }
    ///
    /// Check the Dbg struct is valid.
    ///
    /// * no parallel edge
    /// * no duplicated node
    /// * copy_nums of nodes are consistent
    /// * copy_nums of edges are consistent
    /// * backend DiGraph is valid
    ///
    pub fn is_valid(&self) -> bool {
        self.has_consistent_node_copy_nums()
            && self.has_consistent_edge_copy_nums()
            && self.has_no_duplicated_node()
            && self.has_no_parallel_edge()
            && self.is_graph_valid()
    }
}

///
/// For Node/Edge CopyNums vector
///
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    /// Create the NodeVec with copy numbers of each node
    ///
    pub fn to_node_copy_nums(&self) -> NodeCopyNums {
        let mut v: NodeCopyNums = NodeCopyNums::new(self.n_nodes(), 0);
        for (node, weight) in self.nodes() {
            v[node] = weight.copy_num();
        }
        v
    }
    ///
    /// Get a vector of edge copy numbers (`vec[edge.index()] = edge.copy_num()`)
    ///
    /// copy_num of warp edge (tail -> head) is unassigned, so it will be
    /// filled with 0.
    ///
    pub fn to_edge_copy_nums(&self) -> Option<EdgeCopyNums> {
        if self.is_edge_copy_nums_assigned() {
            let mut v: EdgeCopyNums = EdgeCopyNums::new(self.n_edges(), 0);
            for (edge, _, _, weight) in self.edges() {
                if self.is_warp_edge(edge) {
                    v[edge] = 0;
                } else {
                    v[edge] = weight.copy_num().unwrap();
                }
            }
            Some(v)
        } else {
            None
        }
    }
    ///
    /// Assign copy numbers to all nodes at a time, specified by copy_nums NodeCopyNums vector.
    ///
    /// It returns the copy_nums are updated or not.
    ///
    pub fn set_node_copy_nums(&mut self, copy_nums: &NodeCopyNums) -> bool {
        // vector length assertion
        assert!(copy_nums.len() == self.n_nodes());
        let mut is_updated = false;

        for (i, node_weight_mut) in self.graph.node_weights_mut().enumerate() {
            let node = NodeIndex::new(i);
            let copy_num = copy_nums[node];
            if node_weight_mut.copy_num() != copy_num {
                is_updated = true;
            }
            node_weight_mut.set_copy_num(copy_num);
        }

        is_updated
    }
    ///
    /// Assign copy numbers to all edges at a time, specified by copy_nums EdgeCopyNums vector.
    ///
    /// It returns the copy_nums are updated or not.
    ///
    pub fn set_edge_copy_nums(&mut self, copy_nums: Option<&EdgeCopyNums>) -> bool {
        // vector length assertion
        assert!(copy_nums.is_none() || copy_nums.unwrap().len() == self.n_edges());
        let mut is_updated = false;

        for edge in self.graph.edge_indices() {
            let is_warp_edge = self.is_warp_edge(edge);
            let edge_weight_mut = self.graph.edge_weight_mut(edge).unwrap();

            let copy_num = if is_warp_edge {
                // copy num of warp edges is always unassigned.
                None
            } else {
                // assign the value of vector for normal edges.
                match copy_nums {
                    None => None,
                    Some(copy_nums) => Some(copy_nums[edge]),
                }
            };

            // check the new_copy_num is different from the original copy num
            if edge_weight_mut.copy_num() != copy_num {
                is_updated = true;
            }
            edge_weight_mut.set_copy_num(copy_num)
        }

        is_updated
    }
    ///
    /// Calculate a path (= a list of nodes) that gives the sequence as a emission along the path
    /// as a sequence with SeqStyle::Linear.
    ///
    fn to_nodes_of_seq(&self, seq: &[u8]) -> Option<Vec<NodeIndex>> {
        let styled_sequence = StyledSequence::new(seq.to_vec(), SeqStyle::Linear);
        self.to_nodes_of_styled_seq(&styled_sequence)
    }
    ///
    /// Given a StyledSequence, calculate a path (= a list of nodes) that gives
    /// the sequence as a emission along the path.
    ///
    /// If the sequence cannot be emitted using a path in the dbg, returns None.
    ///
    fn to_nodes_of_styled_seq(&self, seq: &StyledSequence) -> Option<Vec<NodeIndex>> {
        let m = self.to_kmer_map();
        let mut nodes: Vec<NodeIndex> = Vec::new();
        for kmer in styled_sequence_to_kmers(seq, self.k()) {
            match m.get(&kmer) {
                None => return None,
                Some(&node) => nodes.push(node),
            }
        }
        Some(nodes)
    }
    ///
    /// generate node/edge copy numbers of the given sequence as SeqStyle::Linear
    ///
    pub fn to_copy_nums_of_seq(&self, seq: &[u8]) -> Option<(NodeCopyNums, EdgeCopyNums)> {
        let styled_sequence = StyledSequence::new(seq.to_vec(), SeqStyle::Linear);
        self.to_copy_nums_of_styled_seq(&styled_sequence)
    }
    ///
    /// generate node/edge copy numbers of the given sequence
    ///
    pub fn to_copy_nums_of_styled_seq(
        &self,
        seq: &StyledSequence,
    ) -> Option<(NodeCopyNums, EdgeCopyNums)> {
        // vectors to be returned
        let mut nc: NodeCopyNums = NodeCopyNums::new(self.n_nodes(), 0);
        let mut ec: EdgeCopyNums = EdgeCopyNums::new(self.n_edges(), 0);

        match self.to_nodes_of_styled_seq(seq) {
            None => None,
            Some(nodes) => {
                // add node counts
                for &node in nodes.iter() {
                    nc[node] += 1;
                }

                // add edge counts
                for (&node_a, &node_b) in nodes.iter().tuple_windows() {
                    let edge = self
                        .find_edge(node_a, node_b)
                        .expect("there is no corresponding edge in the dbg");
                    ec[edge] += 1;
                }

                if !seq.style().is_fragment() {
                    // add edge between tail and head if the path is circle.
                    let head = nodes.first().unwrap();
                    let tail = nodes.last().unwrap();
                    let edge = self
                        .find_edge(*tail, *head)
                        .expect("there is no corresponding edge in the dbg");
                    ec[edge] += 1;
                }

                Some((nc, ec))
            }
        }
    }
}

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    /// to kmer count profile
    ///
    pub fn to_kmer_profile(&self) -> HashMap<N::Kmer, CopyNum> {
        let mut hm = HashMap::default();
        for (_node, weight) in self.nodes() {
            hm.insert(weight.kmer().clone(), weight.copy_num());
        }
        hm
    }
}

///
/// Seq addition
/// TODO
///
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    /// find the kmer in the de bruijn graph
    ///
    pub fn get_kmer(&self, kmer: &N::Kmer) -> Option<NodeIndex> {
        self.nodes()
            .find(|(_, weight)| weight.kmer() == kmer)
            .map(|(node, _)| node)
    }
    ///
    /// create hashmap from kmer to node index in the dbg.
    ///
    pub fn to_kmer_map(&self) -> HashMap<N::Kmer, NodeIndex> {
        let mut hm = HashMap::default();
        for (node, weight) in self.nodes() {
            hm.insert(weight.kmer().clone(), node);
        }
        hm
    }
    ///
    /// WIP
    ///
    /// add kmer to the de bruijn graph, if not exists.
    ///
    /// # TODOs
    ///
    /// * determine the correct behaviour when the same kmer exists?
    /// * after this addition, the copy-number consistency will be broken.
    ///
    pub fn add_kmer(&mut self, kmer: N::Kmer, copy_num: CopyNum) -> Option<NodeIndex> {
        match self.get_kmer(&kmer) {
            Some(node) => {
                panic!("kmer {} is already exists as node {:?}", kmer, node);
            }
            None => {
                let node = self.graph.add_node(N::new(kmer, copy_num));
                // TODO add edges between parents/childs
                Some(node)
            }
        }
    }
}

///
/// Basic constructors
///
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    /// plain constructor of dbg
    pub fn from_digraph(k: usize, graph: DiGraph<N, E>) -> Self {
        Dbg { k, graph }
    }
    /// Create an empty de bruijn graph with no nodes and edges.
    pub fn empty(k: usize) -> Self {
        Dbg {
            k,
            graph: DiGraph::new(),
        }
    }
    /// Convert HashDbg<K> into Dbg
    pub fn from_hashdbg(d: &HashDbg<N::Kmer>) -> Self {
        let mut graph = DiGraph::new();
        // a temporary map from Kmer to NodeIndex
        let mut ids: HashMap<N::Kmer, NodeIndex> = HashMap::default();

        // (1) add a node for each kmer
        for kmer in d.kmers() {
            let node = graph.add_node(N::new(kmer.clone(), d.get(kmer)));
            ids.insert(kmer.clone(), node);
        }

        // (2) add an edge from kmer to its childs
        for kmer in d.kmers() {
            let v = *ids.get(kmer).unwrap();
            for child in d.childs(kmer) {
                let w = *ids.get(&child).unwrap();
                graph.add_edge(v, w, E::new(None));
            }
        }

        Self::from_digraph(d.k(), graph)
    }
    /// Construct Dbg from a Sequence via converting HashDbg into Dbg.
    pub fn from_seq(k: usize, seq: &[u8]) -> Self {
        let hd = HashDbg::from_seq(k, seq);
        Self::from_hashdbg(&hd)
    }
    /// Construct Dbg from Reads
    /// via converting HashDbg into Dbg.
    pub fn from_reads(k: usize, reads: &Reads) -> Self {
        let hd = HashDbg::from_reads(k, reads);
        Self::from_hashdbg(&hd)
    }
    ///
    /// Convert into edge-centric de bruijn graph
    ///
    pub fn to_edbg(&self) -> SimpleEDbg<N::Kmer> {
        self.to_edbg_with_attr(None)
    }
    ///
    /// Convert into edge-centric de bruijn graph with attributes
    ///
    /// if no attributes vector is given, the default value of the type
    /// will be assigned to `edge.attribute`.
    ///
    /// EdgeIndex of edges in edbg corresponds to NodeIndex of nodes in dbg.
    ///
    pub fn to_edbg_with_attr<A: Copy + PartialEq + Default>(
        &self,
        attrs: Option<&NodeVec<DenseStorage<A>>>,
    ) -> SimpleEDbgWithAttr<N::Kmer, A> {
        let mut graph = DiGraph::new();
        let mut nodes: HashMap<N::Kmer, NodeIndex> = HashMap::default();

        for (node, weight) in self.nodes() {
            let kmer = weight.kmer().clone();
            let copy_num = weight.copy_num();

            // add prefix node if not exists
            let prefix = kmer.prefix();
            let v = match nodes.get(&prefix) {
                None => {
                    let node = graph.add_node(SimpleEDbgNode::new(prefix.clone()));
                    nodes.insert(prefix, node);
                    node
                }
                Some(&node) => node,
            };

            // add suffix node if not exists
            let suffix = kmer.suffix();
            let w = match nodes.get(&suffix) {
                None => {
                    let node = graph.add_node(SimpleEDbgNode::new(suffix.clone()));
                    nodes.insert(suffix, node);
                    node
                }
                Some(&node) => node,
            };

            // add an edge for this kmer
            let attr = match attrs {
                Some(attrs) => attrs[node],
                None => A::default(),
            };
            graph.add_edge(
                v,
                w,
                SimpleEDbgEdgeWithAttr::new_with_attr(kmer, copy_num, node, attr),
            );
        }
        SimpleEDbgWithAttr::new(self.k(), graph)
    }
    ///
    /// Create a `k+1` dbg from the `k` dbg whose edge copy numbers are
    /// consistently assigned.
    ///
    /// * Nodes in k+1-dbg
    ///     an edge in k-dbg represents a k+1-mer, and it will correspond to a node in k+1-dbg.
    /// * Edges in k+1-dbg
    ///     a pair of edges (e1, e2) whose end-point and start-point is the same, has an edge in
    ///     k+1-dbg.
    ///
    /// ## TODOs
    ///
    /// * NNNN should be treated specially. (required)
    /// * zero-copy-number nodes should be purged. (not required)
    ///
    pub fn to_kp1_dbg(&self) -> Dbg<N, E> {
        assert!(self.is_edge_copy_nums_assigned());
        assert!(self.has_consistent_edge_copy_nums());

        let mut graph = DiGraph::new();

        // mapping from "edge in k-dbg" into "node in k+1-dbg".
        let mut ids: HashMap<EdgeIndex, NodeIndex> = HashMap::default();

        // a edge in k-dbg is corresponds to a node in k+1-dbg.
        for (edge, s, t, weight) in self.edges() {
            if !self.is_warp_edge(edge) {
                let kmer = self.kmer(s).join(self.kmer(t));
                let copy_num = weight.copy_num().unwrap_or_else(|| {
                    panic!(
                        "edge {} have no copy nums even though it is not warp edge.",
                        edge.index()
                    )
                });
                let node = graph.add_node(N::new(kmer, copy_num));
                ids.insert(edge, node);
            }
        }

        for (node, weight) in self.nodes() {
            match (weight.is_head(), weight.is_tail()) {
                (false, false) => {
                    // add an edge between all in-edges and out-edges pair.
                    // --e1--> node --e2--> in k-dbg
                    for (e1, _, _) in self.parents(node) {
                        let v1 = ids.get(&e1).unwrap();
                        for (e2, _, _) in self.childs(node) {
                            let v2 = ids.get(&e2).unwrap();
                            // copy numbers of edges in k+1 is ambiguous.
                            graph.add_edge(*v1, *v2, E::new(None));
                        }
                    }
                }
                _ => {}
            }
        }

        // intersections
        let tips = self.tips();
        let in_nodes: Vec<NodeIndex> = tips
            .iter_in_nodes()
            .map(|v| {
                // add a node of in_node
                graph.add_node(N::new(
                    self.node(v).kmer().extend_tail(),
                    self.node(v).copy_num(),
                ))
            })
            .collect();
        let out_nodes: Vec<NodeIndex> = tips
            .iter_out_nodes()
            .map(|v| {
                // add a node of out_node
                graph.add_node(N::new(
                    self.node(v).kmer().extend_head(),
                    self.node(v).copy_num(),
                ))
            })
            .collect();
        for (&w1, &w2) in iproduct!(in_nodes.iter(), out_nodes.iter()) {
            graph.add_edge(w1, w2, E::new(None));
        }

        // add an edge for a in_edges of tail nodes
        for (v, &w) in izip!(tips.iter_in_nodes(), in_nodes.iter()) {
            for (e, _, _) in self.parents(v) {
                let w2 = ids.get(&e).unwrap();
                // w2: parent(YXNN) -> w: tail(XNNN)
                graph.add_edge(*w2, w, E::new(None));
            }
        }

        // add an edge for a out_edges of head nodes
        for (v, &w) in izip!(tips.iter_out_nodes(), out_nodes.iter()) {
            for (e, _, _) in self.childs(v) {
                let w2 = ids.get(&e).unwrap();
                // w: head(NNNX) -> w2: child(NNXY)
                graph.add_edge(w, *w2, E::new(None));
            }
        }

        Self::from_digraph(self.k() + 1, graph)
    }
}

#[cfg(test)]
mod tests {
    use super::super::mocks::*;
    use super::*;
    use crate::common::ni;
    use crate::common::sequence_to_string;
    use crate::dbg::edge_centric::EDbgEdge;
    use crate::kmer::veckmer::VecKmer;

    #[test]
    fn dbg_new() {
        let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"ATCGGCT");
        println!("{}", hd);
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_hashdbg(&hd);
        println!("{}", dbg);
        for (node, weight) in dbg.nodes() {
            println!("{:?} {}", node, weight);
            for (edge, child, weight) in dbg.childs(node) {
                println!("{:?} {:?} {}", edge, child, weight);
            }
            for (edge, parent, weight) in dbg.parents(node) {
                println!("{:?} {:?} {}", edge, parent, weight);
            }
        }
    }
    #[test]
    fn dbg_to_edbg() {
        let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"ATCGGCT");
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_hashdbg(&hd);
        println!("{}", dbg);
        let edbg = dbg.to_edbg();
        println!("{}", edbg);
        assert_eq!(edbg.n_edges(), dbg.n_nodes());
    }
    #[test]
    fn dbg_to_edbg_with_attr() {
        let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"ATCGGCT");
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_hashdbg(&hd);

        let a = 10.1;
        let b = 1.1;
        let mut v: NodeVec<DenseStorage<f64>> = NodeVec::new(dbg.n_nodes(), a);
        v[ni(1)] = b;

        let edbg = dbg.to_edbg_with_attr(Some(&v));
        println!("{}", edbg);
        for (edge, v, w, weight) in edbg.edges() {
            if weight.origin_node() == ni(1) {
                assert_eq!(weight.attribute, b);
            } else {
                assert_eq!(weight.attribute, a);
            }
        }
    }
    #[test]
    fn dbg_kmer() {
        let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"ATCGGCT");
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_hashdbg(&hd);
        println!("{}", dbg);
        assert!(!dbg.is_edge_copy_nums_assigned());
        println!("{:?}", dbg.get_kmer(&VecKmer::from_bases(b"ATCG")));
        assert_eq!(dbg.get_kmer(&VecKmer::from_bases(b"ATCA")), None);
        assert_eq!(dbg.get_kmer(&VecKmer::from_bases(b"ATCG")), Some(ni(9)));

        let m = dbg.to_kmer_map();
        assert_eq!(m.get(&VecKmer::from_bases(b"ATCA")).copied(), None);
        assert_eq!(m.get(&VecKmer::from_bases(b"TCGG")).copied(), Some(ni(5)));
        println!("{:?}", m);
    }
    #[test]
    fn dbg_copy_numbers() {
        let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"ATCGGCT");
        let mut dbg: SimpleDbg<VecKmer> = SimpleDbg::from_hashdbg(&hd);
        println!("{}", dbg);
        assert_eq!(dbg.to_node_copy_nums().to_vec(), vec![1; dbg.n_nodes()]);
        assert_eq!(dbg.to_edge_copy_nums(), None);
        assert_eq!(dbg.is_edge_copy_nums_assigned(), false);

        // node copy numbers assignment
        dbg.set_node_copy_nums(&NodeCopyNums::from_slice(&vec![0; dbg.n_nodes()], 0));
        assert_eq!(dbg.to_node_copy_nums().to_vec(), vec![0; dbg.n_nodes()],);

        // edge copy numbers assignment
        dbg.set_edge_copy_nums(Some(&EdgeCopyNums::from_slice(&vec![1; dbg.n_edges()], 0)));
        assert_eq!(
            dbg.to_edge_copy_nums().unwrap().to_vec(),
            vec![1, 1, 1, 1, 1, 1, 0, 1, 1, 1]
        );
        assert_eq!(dbg.is_edge_copy_nums_assigned(), true);

        let nodes = dbg.to_nodes_of_seq(b"ATCGGCT").unwrap();
        println!("nodes={:?}", nodes);
        assert_eq!(
            nodes,
            vec![
                ni(7),
                ni(3),
                ni(4),
                ni(9),
                ni(5),
                ni(8),
                ni(0),
                ni(1),
                ni(2),
                ni(6)
            ]
        );
        println!("{}", sequence_to_string(&dbg.path_as_sequence(&nodes)));
        assert_eq!(dbg.path_as_sequence(&nodes), b"ATCGGCTNNN");

        let (ncn, ecn) = dbg.to_copy_nums_of_seq(b"ATCGGCT").unwrap();
        assert_eq!(ncn.to_vec(), vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);
        assert_eq!(ecn.to_vec(), vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);

        assert_eq!(dbg.has_consistent_node_copy_nums(), true);
        assert_eq!(dbg.has_consistent_edge_copy_nums(), true);

        dbg.set_edge_copy_nums(Some(&EdgeCopyNums::from_slice(
            &vec![1, 1, 1, 3, 1, 1, 3, 1, 1, 1],
            0,
        )));
        assert_eq!(dbg.has_consistent_edge_copy_nums(), false);

        dbg.set_edge_copy_nums(Some(&EdgeCopyNums::from_slice(
            &vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            0,
        )));
        assert_eq!(dbg.has_consistent_edge_copy_nums(), true);

        dbg.set_node_copy_nums(&NodeCopyNums::from_slice(
            &vec![1, 1, 1, 3, 1, 1, 3, 1, 1, 1],
            0,
        ));
        assert_eq!(dbg.has_consistent_node_copy_nums(), false);
    }
    #[test]
    fn manual_dbg() {
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::empty(4);
        // dbg.add_node();
        // dbg.add_seq(b"ATCGAT");
    }
    #[test]
    fn dbg_tips() {
        // base
        let dbg = mock_base();
        let tips = dbg.tips();
        println!("{}", dbg);
        println!("{}", tips);
        assert_eq!(tips.n_in_nodes(), 1);
        assert_eq!(tips.in_node(0), ni(6));
        assert_eq!(tips.n_out_nodes(), 1);
        assert_eq!(tips.out_node(0), ni(7));
        assert!(tips.km1mer().is_null());

        // intersection
        let dbg = mock_intersection();
        let tips = dbg.tips();
        println!("{}", dbg);
        println!("{}", tips);
        assert_eq!(tips.n_in_nodes(), 2);
        assert_eq!(tips.in_node(0), ni(3));
        assert_eq!(tips.in_node(1), ni(7));
        assert_eq!(tips.n_out_nodes(), 2);
        assert_eq!(tips.out_node(0), ni(8));
        assert_eq!(tips.out_node(1), ni(14));
        assert!(tips.km1mer().is_null());
    }
    #[test]
    fn dbg_extension() {
        let mut dbg = mock_intersection();
        // println!("{}", dbg);
        let (ncn, ecn) = dbg.to_copy_nums_of_seq(b"AACTAGCTT").unwrap();
        dbg.set_node_copy_nums(&ncn);
        dbg.set_edge_copy_nums(Some(&ecn));
        println!("{}", dbg);
        // println!("{}", dbg.to_cytoscape());
        let dbg2 = dbg.to_kp1_dbg();
        println!("{}", dbg2);
        assert!(dbg2.has_consistent_node_copy_nums());

        let seqs = dbg.to_seqs();
        for seq in seqs.iter() {
            println!("dbg={}", sequence_to_string(seq));
        }

        let seqs = dbg2.to_seqs();
        for seq in seqs.iter() {
            println!("dbg2={}", sequence_to_string(seq));
        }
    }
    #[test]
    fn dbg_extension_intersection_a() {
        //
        // [pattern a]
        // representing ATAGCT and TAAGCC
        //
        let mut dbg = mock_intersection_small();
        println!("{}", dbg);
        println!("{}", dbg.to_dot());
        assert!(dbg.is_valid());
        assert!(!dbg.is_edge_copy_nums_assigned());
        assert_eq!(format!("{}", dbg), "4,L:TAAGCT,L:ATAGCC");

        // discarded edges
        // * TAGC 17 -> AGCC 5
        // * AAGC 11 -> AGCT 12
        // * warp edges 8->9, 8->3, 6->9, 6->3
        let ecn = EdgeCopyNums::from_slice(
            &[
                1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1,
            ],
            0,
        );
        let mut ecn2 = EdgeCopyNums::new(dbg.n_edges(), 1);
        ecn2[dbg
            .find_edge_from_kp1mer(&VecKmer::from_bases(b"TAGCC"))
            .unwrap()] = 0;
        ecn2[dbg
            .find_edge_from_kp1mer(&VecKmer::from_bases(b"AAGCT"))
            .unwrap()] = 0;
        println!("{}", ecn);
        println!("{}", ecn2);
        dbg.set_edge_copy_nums(Some(&ecn));
        assert!(dbg.is_edge_copy_nums_assigned());
        assert_eq!(dbg.to_edge_copy_nums().unwrap(), ecn);

        let dbg2 = dbg.to_kp1_dbg();
        println!("{}", dbg2);
        println!("{}", dbg2.to_dot());
        assert!(dbg2.is_valid());
        assert!(!dbg2.is_edge_copy_nums_assigned());
        assert_eq!(format!("{}", dbg2), "5,L:TAAGCC,L:ATAGCT");
    }
    #[test]
    fn dbg_extension_intersection_b() {
        //
        // [pattern b]
        // representing ATAGCC and TAAGCT
        //
        let mut dbg = mock_intersection_small();
        println!("{}", dbg);
        println!("{}", dbg.to_dot());
        assert!(dbg.is_valid());
        assert!(!dbg.is_edge_copy_nums_assigned());
        assert_eq!(format!("{}", dbg), "4,L:TAAGCT,L:ATAGCC");

        // discarded edges
        // * TAGC 17 -> AGCT 12
        // * AAGC 11 -> AGCC 5
        // * warp edges 8->9, 8->3, 6->9, 6->3
        let ecn = EdgeCopyNums::from_slice(
            &[
                1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0,
            ],
            0,
        );
        println!("{}", ecn);
        dbg.set_edge_copy_nums(Some(&ecn));
        assert!(dbg.is_edge_copy_nums_assigned());
        assert_eq!(dbg.to_edge_copy_nums().unwrap(), ecn);

        let dbg2 = dbg.to_kp1_dbg();
        println!("{}", dbg2);
        println!("{}", dbg2.to_dot());
        assert!(dbg2.is_valid());
        assert!(!dbg2.is_edge_copy_nums_assigned());
        assert_eq!(format!("{}", dbg2), "5,L:TAAGCT,L:ATAGCC");
    }
    #[test]
    fn dbg_extension_2() {
        let mut dbg = mock_manual();
        println!("{}", dbg);
        println!("{}", dbg.to_dot());
        assert!(dbg.is_valid());
        assert!(!dbg.is_edge_copy_nums_assigned());
        assert_eq!(format!("{}", dbg), "4,L:ATCGAGCATG");

        // (1) set the edge_copy_nums
        let mut ecn: EdgeCopyNums =
            EdgeCopyNums::from_slice(&[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0], 0);
        dbg.set_edge_copy_nums(Some(&ecn));
        assert!(dbg.is_edge_copy_nums_assigned());

        // converting back gives the same vector?
        let ecn2 = dbg.to_edge_copy_nums().unwrap();
        println!("{}", ecn2);
        assert_eq!(ecn, ecn2);

        // the edges have proper copy nums?
        for (edge, _, _, weight) in dbg.edges() {
            if dbg.is_warp_edge(edge) {
                assert!(weight.copy_num().is_none());
            } else {
                assert!(weight.copy_num().is_some());
                assert_eq!(weight.copy_num().unwrap(), 1);
            }
        }

        // (2) upgrade the dbg
        let dbg2 = dbg.to_kp1_dbg();
        println!("{}", dbg2);
        println!("{}", dbg2.to_dot());
        assert!(dbg2.is_valid());
        assert!(!dbg2.is_edge_copy_nums_assigned());

        assert_eq!(format!("{}", dbg2), "5,L:ATCGAGCATG");
    }
    #[test]
    fn dbg_clone() {
        let mut g: DiGraph<u32, ()> = DiGraph::new();
        let v = g.add_node(10);
        println!("{}", g.node_count());
        println!("{}", g.edge_count());

        let g2 = (&g).clone();
        let v = g.add_node(11);

        println!("{}", g.node_count());
        println!("{}", g.edge_count());

        println!("{}", g2.node_count());
        println!("{}", g2.edge_count());
    }
}
