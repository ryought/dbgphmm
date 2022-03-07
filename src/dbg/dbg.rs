//!
//! De bruijn graph definitions
//!
//!
use super::edge_centric::{SimpleEDbg, SimpleEDbgEdge, SimpleEDbgNode};
use super::impls::{SimpleDbg, SimpleDbgEdge, SimpleDbgNode};
use crate::common::{CopyNum, Sequence};
use crate::dbg::hashdbg_v2::HashDbg;
use crate::graph::iterators::{ChildEdges, EdgesIterator, NodesIterator, ParentEdges};
use crate::kmer::kmer::{sequence_to_kmers, Kmer, KmerLike};
use crate::vector::{DenseStorage, EdgeVec, NodeVec};
use fnv::FnvHashMap as HashMap;
use itertools::Itertools;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

pub type NodeCopyNums = NodeVec<DenseStorage<CopyNum>>;
pub type EdgeCopyNums = EdgeVec<DenseStorage<CopyNum>>;

///
/// (Node-centric) De bruijn graph struct
/// k
///
pub struct Dbg<N: DbgNode, E: DbgEdge> {
    k: usize,
    pub graph: DiGraph<N, E>,
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
    /// Modify copy number count of this node
    fn set_copy_num(&mut self, copy_num: CopyNum);
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
}

///
/// Basic properties
///
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    /// CopyNums of nodes are consistent, that is
    /// 'sum of copynums of childs' == 'sum of copynums of siblings'
    /// for each kmers
    fn has_consistent_node_copy_nums(&self) -> bool {
        unimplemented!();
    }
    fn has_consistent_edge_copy_nums(&self) -> bool {
        unimplemented!();
    }
    /// Check if all the edges has a copy number
    fn is_edge_copy_nums_assigned(&self) -> bool {
        self.edges().all(|(_, _, _, ew)| ew.copy_num().is_some())
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
        // TODO assert that node copy nums are consistent
        let mut v: NodeCopyNums = NodeCopyNums::new(self.n_nodes(), 0);
        for (node, weight) in self.nodes() {
            v[node] = weight.copy_num();
        }
        v
    }
    ///
    /// Get a vector of edge copy numbers (`vec[edge.index()] = edge.copy_num()`)
    ///
    pub fn to_edge_copy_nums(&self) -> Option<EdgeCopyNums> {
        // TODO assert edge copy nums are consistent
        if self.is_edge_copy_nums_assigned() {
            let mut v: EdgeCopyNums = EdgeCopyNums::new(self.n_edges(), 0);
            for (edge, _, _, weight) in self.edges() {
                v[edge] = weight.copy_num().unwrap();
            }
            Some(v)
        } else {
            None
        }
    }
    ///
    /// Assign copy numbers to all nodes at a time, specified by copy_nums NodeCopyNums vector.
    ///
    pub fn set_node_copy_nums(&mut self, copy_nums: &NodeCopyNums) {
        for (i, node_weight_mut) in self.graph.node_weights_mut().enumerate() {
            let node = NodeIndex::new(i);
            node_weight_mut.set_copy_num(copy_nums[node])
        }
    }
    ///
    /// Assign copy numbers to all edges at a time, specified by copy_nums EdgeCopyNums vector.
    ///
    pub fn set_edge_copy_nums(&mut self, copy_nums: Option<&EdgeCopyNums>) {
        for (i, edge_weight_mut) in self.graph.edge_weights_mut().enumerate() {
            let edge = EdgeIndex::new(i);
            let copy_num = match copy_nums {
                None => None,
                Some(copy_nums) => Some(copy_nums[edge]),
            };
            edge_weight_mut.set_copy_num(copy_num)
        }
    }
    fn to_nodes_of_seq(&self, seq: &[u8]) -> Option<Vec<NodeIndex>> {
        let m = self.to_kmer_map();
        let mut nodes: Vec<NodeIndex> = Vec::new();
        for kmer in sequence_to_kmers(seq, self.k()) {
            match m.get(&kmer) {
                None => return None,
                Some(&node) => nodes.push(node),
            }
        }
        Some(nodes)
    }
    ///
    /// generate node/edge copy numbers of the given sequence
    ///
    pub fn to_copy_nums_of_seq(&self, seq: &[u8]) -> Option<(NodeCopyNums, EdgeCopyNums)> {
        // vectors to be returned
        let mut nc: NodeCopyNums = NodeCopyNums::new(self.n_nodes(), 0);
        let mut ec: EdgeCopyNums = EdgeCopyNums::new(self.n_edges(), 0);

        match self.to_nodes_of_seq(seq) {
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

                Some((nc, ec))
            }
        }
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
    ///
    /// Convert into edge-centric de bruijn graph
    ///
    pub fn to_edbg(&self) -> SimpleEDbg<N::Kmer> {
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
            graph.add_edge(v, w, SimpleEDbgEdge::new(kmer, copy_num, node));
        }
        SimpleEDbg::new(self.k(), graph)
    }
    ///
    /// Create a `k+1` dbg from the `k` dbg whose edge copy numbers are
    /// consistently assigned.
    ///
    pub fn to_kp1_dbg(&self) -> Dbg<N, E> {
        unimplemented!();
    }
}

#[cfg(test)]
mod tests {
    use super::super::mocks::*;
    use super::*;
    use crate::common::ni;
    use crate::common::sequence_to_string;
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

        // node copy numbers assignment
        dbg.set_node_copy_nums(&NodeCopyNums::from_slice(&vec![0; dbg.n_nodes()], 0));
        assert_eq!(dbg.to_node_copy_nums().to_vec(), vec![0; dbg.n_nodes()],);

        // edge copy numbers assignment
        dbg.set_edge_copy_nums(Some(&EdgeCopyNums::from_slice(&vec![1; dbg.n_edges()], 0)));
        assert_eq!(
            dbg.to_edge_copy_nums().unwrap().to_vec(),
            vec![1; dbg.n_edges()]
        );

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
    }
    #[test]
    fn manual_dbg() {
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::empty(4);
        // dbg.add_node();
        // dbg.add_seq(b"ATCGAT");
    }
    #[test]
    fn dbg_extension() {
        let dbg = mock_intersection();
        println!("{}", dbg);
    }
}
