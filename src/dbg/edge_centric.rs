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
use crate::kmer::kmer::{Kmer, KmerLike};
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
pub mod impls;
pub use impls::{SimpleEDbg, SimpleEDbgEdge, SimpleEDbgNode};
use petgraph::dot::Dot;

///
/// (Edge-centric) De bruijn graph struct
///
pub struct EDbg<N: EDbgNode, E: EDbgEdge> {
    k: usize,
    pub graph: DiGraph<N, E>,
}

///
/// Trait for nodes in edge-centric dbg `EDbg`
///
pub trait EDbgNode {}

///
/// Trait for edges in edge-centric dbg `EDbg`
///
pub trait EDbgEdge {
    type Kmer: KmerLike;
    ///
    /// k-mer of this edge of the EDbg
    fn kmer(&self) -> &Self::Kmer;
    ///
    /// Copy number count of this edge in EDbg
    fn copy_num(&self) -> CopyNum;
}

///
/// Basic graph operations for Dbg
///
impl<N: EDbgNode, E: EDbgEdge> EDbg<N, E> {
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
    /// plain constructor of edbg
    pub fn new(k: usize, graph: DiGraph<N, E>) -> Self {
        EDbg { k, graph }
    }
}

impl<N, E> std::fmt::Display for EDbg<N, E>
where
    N: EDbgNode + std::fmt::Display,
    E: EDbgEdge + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.graph, &[]))
    }
}
