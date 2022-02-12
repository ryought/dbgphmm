//!
//! De bruijn graph definitions
//!
//!
use super::impls::{SimpleDbg, SimpleDbgEdge, SimpleDbgNode};
use crate::common::{CopyNum, Sequence};
use crate::dbg::hashdbg_v2::HashDbg;
use crate::graph::iterators::{ChildEdges, EdgesIterator, NodesIterator, ParentEdges};
use crate::kmer::kmer::{Kmer, KmerLike};
use crate::vector::{DenseStorage, EdgeVec, NodeVec};
use fnv::FnvHashMap as HashMap;
use petgraph::dot::Dot;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

pub type NodeCopyNums = NodeVec<DenseStorage<CopyNum>>;
pub type EdgeCopyNums = EdgeVec<DenseStorage<CopyNum>>;

///
/// (Node-centric) De bruijn graph struct
/// k
///
pub struct Dbg<N: DbgNode, E: DbgEdge> {
    k: usize,
    graph: DiGraph<N, E>,
}

///
/// Trait for nodes in Dbg
///
pub trait DbgNode {
    type Kmer: KmerLike;
    fn new(kmer: Self::Kmer, copy_num: CopyNum) -> Self;
    ///
    /// Kmer of this node of the Dbg
    fn kmer(&self) -> &Self::Kmer;
    ///
    /// Copy number count of this node in Dbg
    fn copy_num(&self) -> CopyNum;
    ///
    /// Single base assigned to this node in Dbg
    /// Last base of kmer will be used as an emission
    fn emission(&self) -> u8 {
        self.kmer().last()
    }
}

///
/// Trait for edges in Dbg
///
pub trait DbgEdge {
    fn new(copy_num: Option<CopyNum>) -> Self;
    fn copy_num(&self) -> Option<CopyNum>;
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
    /// Return the number of nodes in the graph
    pub fn n_nodes(&self) -> usize {
        self.graph.node_count()
    }
    /// Return the number of edges in the graph
    pub fn n_edges(&self) -> usize {
        self.graph.edge_count()
    }
    /// convert node to emission
    pub fn emission(&self, node: NodeIndex) -> u8 {
        self.graph.node_weight(node).unwrap().emission()
    }
    /// check if two nodes `a, b: NodeIndex` is connected or not
    pub fn contains_edge(&self, a: NodeIndex, b: NodeIndex) -> bool {
        self.graph.contains_edge(a, b)
    }
}

///
/// Basic properties
///
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    pub fn is_consistent(&self) {
        unimplemented!();
    }
}

///
/// Basic constructors
///
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    /// plain constructor of dbg
    pub fn new(k: usize, graph: DiGraph<N, E>) -> Self {
        Dbg { k, graph }
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

        Self::new(d.k(), graph)
    }
}

impl<N, E> std::fmt::Display for Dbg<N, E>
where
    N: DbgNode + std::fmt::Display,
    E: DbgEdge + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.graph, &[]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::veckmer::VecKmer;

    #[test]
    fn dbg_new() {
        let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"AAAGCTTGATT");
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
}
