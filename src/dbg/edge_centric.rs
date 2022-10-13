//!
//! Edge-centric de Bruijn graph implementations
//!
//! * Node = k-1-mer overlap
//! * Edge = k-mer
//!
//! # Features
//!
//! * visualization by cytoscape
//! * min-flow optimization
//! * quick-check of copy number consistency
//!
use crate::common::CopyNum;
use crate::graph::iterators::{ChildEdges, EdgesIterator, NodesIterator, ParentEdges};
use crate::hmmv2::trans_table::EdgeFreqs;
use crate::kmer::kmer::{Kmer, KmerLike};
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph::Direction;
pub mod impls;
pub mod output;
use crate::dbg::flow_intersection::{FlowIntersection, FlowIntersectionEdge, FlowIntersectionNode};
pub use impls::{
    SimpleEDbg, SimpleEDbgEdge, SimpleEDbgEdgeWithAttr, SimpleEDbgNode, SimpleEDbgWithAttr,
};
use itertools::iproduct;

///
/// (Edge-centric) De bruijn graph struct
///
pub struct EDbg<N, E> {
    k: usize,
    pub graph: DiGraph<N, E>,
}

///
/// Trait for nodes in edge-centric dbg `EDbg`
///
pub trait EDbgNode {
    type Kmer: KmerLike;
    ///
    /// k-1-mer overlap of this node
    fn km1mer(&self) -> &Self::Kmer;
}

///
/// Fundamental trait for edges in edge-centric dbg `EDbg`
/// If you need copy numbers, use EDbgEdge
///
pub trait EDbgEdgeBase {
    type Kmer: KmerLike;
    ///
    /// k-mer of this edge of the EDbg
    fn kmer(&self) -> &Self::Kmer;
    ///
    /// Index of the node in (node-centric) dbg which this edge is
    /// originated from.
    fn origin_node(&self) -> NodeIndex;
}

///
/// Trait for edges in edge-centric dbg `EDbg`
///
/// EDbgEdgeBase: kmer() and origin_node()
/// EDbgEdge: copy_num()
///
pub trait EDbgEdge: EDbgEdgeBase {
    ///
    /// Copy number count of this edge in EDbg
    fn copy_num(&self) -> CopyNum;
}

///
/// Basic graph operations for Dbg
///
impl<N, E: EDbgEdgeBase> EDbg<N, E> {
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
    /// Return the number of nodes in the graph
    pub fn n_nodes(&self) -> usize {
        self.graph.node_count()
    }
    /// Return the number of edges in the graph
    pub fn n_edges(&self) -> usize {
        self.graph.edge_count()
    }
    /// check if two nodes `a, b: NodeIndex` is connected or not
    pub fn contains_edge(&self, a: NodeIndex, b: NodeIndex) -> bool {
        self.graph.contains_edge(a, b)
    }
}
impl<N: EDbgNode, E: EDbgEdge> EDbg<N, E> {
    /// convert a node into an intersection information
    pub fn intersection(&self, node: NodeIndex) -> IntersectionBase<N::Kmer> {
        let node_weight = self
            .graph
            .node_weight(node)
            .expect("node is not in the graph");

        // list of in/out node indexes
        let in_nodes: Vec<NodeIndex> = self
            .graph
            .edges_directed(node, Direction::Incoming)
            .map(|e| e.weight().origin_node())
            .collect();
        let out_nodes: Vec<NodeIndex> = self
            .graph
            .edges_directed(node, Direction::Outgoing)
            .map(|e| e.weight().origin_node())
            .collect();

        IntersectionBase::new(node_weight.km1mer().clone(), in_nodes, out_nodes)
    }
}

pub struct IntersectionBase<K: KmerLike> {
    km1mer: K,
    in_nodes: Vec<NodeIndex>,
    out_nodes: Vec<NodeIndex>,
}

impl<K: KmerLike> IntersectionBase<K> {
    pub fn new(km1mer: K, in_nodes: Vec<NodeIndex>, out_nodes: Vec<NodeIndex>) -> Self {
        IntersectionBase {
            km1mer,
            in_nodes,
            out_nodes,
        }
    }
    pub fn in_nodes(&self) -> &[NodeIndex] {
        &self.in_nodes
    }
    pub fn out_nodes(&self) -> &[NodeIndex] {
        &self.out_nodes
    }
    pub fn km1mer(&self) -> &K {
        &self.km1mer
    }
}

impl<K: KmerLike> std::fmt::Display for IntersectionBase<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "IntersectionBase({}) in:{:?} out:{:?}",
            self.km1mer(),
            self.in_nodes(),
            self.out_nodes()
        )
    }
}

///
/// Basic properties
///
impl<N: EDbgNode, E: EDbgEdge> EDbg<N, E> {
    ///
    /// check if the CopyNum of edges are consistent,
    /// i.e., the sum of copy numbers of in-edges and out-edges
    /// are the same (for all nodes).
    ///
    pub fn has_consistent_copy_nums(&self) -> bool {
        self.nodes().all(|(node, _)| {
            let in_copy_nums: CopyNum = self
                .graph
                .edges_directed(node, Direction::Incoming)
                .map(|e| e.weight().copy_num())
                .sum();
            let out_copy_nums: CopyNum = self
                .graph
                .edges_directed(node, Direction::Outgoing)
                .map(|e| e.weight().copy_num())
                .sum();
            in_copy_nums == out_copy_nums
        })
    }
}

impl<N, E> EDbg<N, E> {
    /// plain constructor of edbg
    pub fn new(k: usize, graph: DiGraph<N, E>) -> Self {
        EDbg { k, graph }
    }
}
