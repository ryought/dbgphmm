//!
//! Definition of PHMM
//!
use crate::common::{CopyNum, NULL_BASE};
use crate::graph::active_nodes::ActiveNodes;
use crate::graph::iterators::{
    ActiveNodesIterator, ChildEdges, EdgesIterator, NodesIterator, ParentEdges,
};
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;
use petgraph::dot::Dot;
use petgraph::graph::DiGraph;
pub use petgraph::graph::{EdgeIndex, NodeIndex};
use sparsevec::SparseVec;

//
// traits
//

///
/// attribute of PHMM Node
///
/// * `emission(&self) -> u8`
///     The emission assigned to the HMM node.
///
/// * `is_emittable(&self) -> bool`
///     Check if this node has emission or not.
///
/// * `init_prob(&self) -> Prob`
///     The initial transition probability (i.e. Begin state to this node)
///
pub struct PHMMNode {
    ///
    ///
    init_prob: Prob,
    ///
    ///
    emission: u8,
}

impl PHMMNode {
    ///
    ///
    ///
    pub fn new(init_prob: Prob, emission: u8) -> PHMMNode {
        PHMMNode {
            init_prob,
            emission,
        }
    }
    ///
    ///
    ///
    fn emission(&self) -> u8 {
        self.emission
    }
    ///
    /// This node has a valid emission or not.
    ///
    fn is_emittable(&self) -> bool {
        self.emission() != NULL_BASE
    }
    ///
    /// Initial transition probability from Begin state to the node.
    ///
    fn init_prob(&self) -> Prob {
        self.init_prob
    }
}

pub struct PHMMEdge {
    ///
    ///
    ///
    trans_prob: Prob,
}

///
/// PHMMEdge
///
impl PHMMEdge {
    ///
    ///
    ///
    pub fn new(trans_prob: Prob) -> PHMMEdge {
        PHMMEdge { trans_prob }
    }
    ///
    /// Transition probability from the source node to the target node
    /// of this edge.
    ///
    fn trans_prob(&self) -> Prob {
        self.trans_prob
    }
}

impl std::fmt::Display for PHMMNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.is_emittable() {
            write!(f, "{} (p={})", self.emission as char, self.init_prob)
        } else {
            write!(f, "not_emittable (p={})", self.init_prob)
        }
    }
}

impl std::fmt::Display for PHMMEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "p={}", self.trans_prob)
    }
}

/// Profile HMM model
///
pub struct PHMMModel {
    pub param: PHMMParams,
    pub graph: DiGraph<PHMMNode, PHMMEdge>,
}

impl PHMMModel {
    ///
    /// create iterator of all nodes
    /// Item of the iterator is `(NodeIndex, &N)`.
    ///
    pub fn nodes(&self) -> NodesIterator<PHMMNode> {
        NodesIterator::new(&self.graph)
    }
    ///
    /// create iterator of all nodes
    /// Item of the iterator is
    /// `(EdgeIndex, NodeIndex of source, NodeIndex of target, &E)`.
    ///
    pub fn edges(&self) -> EdgesIterator<PHMMEdge> {
        EdgesIterator::new(&self.graph)
    }
    ///
    /// create iterator of all child edges of the node
    ///
    /// Item of the iterator is `(EdgeIndex, NodeIndex, Edge weight)`
    ///
    /// * `EdgeIndex` is index of edge (node, child)
    /// * `NodeIndex` is index of child
    /// * `EdgeWeight` is of the edge transition
    ///
    pub fn childs(&self, node: NodeIndex) -> ChildEdges<PHMMEdge> {
        ChildEdges::new(&self.graph, node)
    }
    ///
    /// create iterator of all parent edges of the node
    ///
    /// Item of the iterator is `(EdgeIndex, NodeIndex, Edge weight)`
    ///
    /// * `EdgeIndex` is index of edge (parent, node)
    /// * `NodeIndex` is index of parent
    /// * `EdgeWeight` is of the edge transition
    ///
    pub fn parents(&self, node: NodeIndex) -> ParentEdges<PHMMEdge> {
        ParentEdges::new(&self.graph, node)
    }
    ///
    /// Return the number of nodes in the graph
    ///
    pub fn n_nodes(&self) -> usize {
        self.graph.node_count()
    }
    ///
    /// Return the number of edges in the graph
    ///
    pub fn n_edges(&self) -> usize {
        self.graph.edge_count()
    }
    ///
    /// Get a reference of node weight
    ///
    pub fn node(&self, node: NodeIndex) -> &PHMMNode {
        self.graph.node_weight(node).unwrap()
    }
    ///
    /// Get a reference of edge weight
    ///
    pub fn edge(&self, edge: EdgeIndex) -> &PHMMEdge {
        self.graph.edge_weight(edge).unwrap()
    }
    ///
    /// Get an index of the edge from a to b.
    ///
    pub fn find_edge(&self, a: NodeIndex, b: NodeIndex) -> Option<EdgeIndex> {
        self.graph.find_edge(a, b)
    }
    //
    // PHMM Specific methods
    //
    ///
    /// Get a node emission
    ///
    pub fn emission(&self, v: NodeIndex) -> u8 {
        self.graph.node_weight(v).unwrap().emission()
    }
    ///
    /// Check that the node is emittable or not.
    ///
    pub fn is_emittable(&self, v: NodeIndex) -> bool {
        self.graph.node_weight(v).unwrap().is_emittable()
    }
    ///
    /// emission probability of observing the emission from
    /// Match state of node v.
    ///
    pub fn p_match_emit(&self, node: NodeIndex, emission: u8) -> Prob {
        if self.emission(node) == emission {
            self.param.p_match
        } else {
            self.param.p_mismatch
        }
    }
    ///
    /// emission probability of observing the emission from
    /// Ins state of node v.
    /// The ret is always equal to `param.p_random`
    ///
    pub fn p_ins_emit(&self) -> Prob {
        self.param.p_random
    }
}

impl std::fmt::Display for PHMMModel {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.graph, &[]))
    }
}

//
// Visualizers
//
impl PHMMModel {
    /// println wrapper of NodeVec, with attributes of each nodes
    pub fn draw_node_vec<T, const M: usize>(&self, nv: &SparseVec<T, NodeIndex, M>)
    where
        T: std::fmt::Display + Copy + PartialOrd,
    {
        for (node, _w) in self.nodes() {
            println!("{:?}({})\t{}", node, self.emission(node) as char, nv[node]);
        }
    }
    /// println wrapper of EdgeVec, with attributes of each edges
    pub fn draw_edge_vec<T, const M: usize>(&self, ev: &SparseVec<T, EdgeIndex, M>)
    where
        T: std::fmt::Display + Copy + PartialOrd,
    {
        for (edge, k, l, _) in self.edges() {
            println!(
                "{:?}({}->{})\t{}",
                edge,
                self.emission(k) as char,
                self.emission(l) as char,
                ev[edge]
            );
        }
    }
}

//
// Tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};
    use crate::graph::mocks::{mock_crossing, mock_linear};
    use crate::graph::seq_graph::SeqGraph;
    use crate::hmmv2::mocks::mock_linear_phmm;
    use crate::prob::p;

    #[test]
    fn phmmmodels_basic_ops() {
        let phmm = mock_linear_phmm(PHMMParams::high_error());
        for (i, (node, weight)) in phmm.nodes().enumerate() {
            assert_eq!(NodeIndex::new(i), node);
            assert_eq!(weight.copy_num, 1);
        }
        for (i, (edge, source, target, weight)) in phmm.edges().enumerate() {
            assert_eq!(EdgeIndex::new(i), edge);
            assert_eq!(NodeIndex::new(i), source);
            assert_eq!(NodeIndex::new(i + 1), target);
        }
    }
    #[test]
    fn phmmmodels_active_nodes() {
        let phmm = mock_linear_phmm(PHMMParams::high_error());
        let n_all: Vec<usize> = (0..10).collect();

        // (1) default nodes iterator
        let n1: Vec<usize> = phmm.nodes().map(|(v, _)| v.index()).collect();
        for node in phmm.nodes() {
            println!("{:?}", node);
        }
        assert_eq!(n1, n_all);

        // (2) all nodes are active nodes
        let n2: Vec<usize> = phmm
            .active_nodes(&ActiveNodes::All)
            .map(|(v, _)| v.index())
            .collect();
        for node in phmm.active_nodes(&ActiveNodes::All) {
            println!("{:?}", node);
        }
        assert_eq!(n2, n_all);

        // (3) three active nodes
        let n3: Vec<usize> = phmm
            .active_nodes(&ActiveNodes::Only(vec![ni(1), ni(5), ni(3)]))
            .map(|(v, _)| v.index())
            .collect();
        assert_eq!(n3, vec![1, 5, 3]);
        for node in phmm.active_nodes(&ActiveNodes::Only(vec![ni(1), ni(5), ni(3)])) {
            println!("{:?}", node);
        }

        // (4) zero active nodes
        let n4: Vec<usize> = phmm
            .active_nodes(&ActiveNodes::Only(Vec::new()))
            .map(|(v, _)| v.index())
            .collect();
        assert_eq!(n4.len(), 0);
    }
    #[test]
    fn hmm_crossing_edge_on() {
        let m1 = mock_crossing(true);
        println!("{}", m1);
        let g1 = m1.to_seq_graph().to_phmm(PHMMParams::default());
        println!("{}", g1);
        assert_eq!(g1.n_nodes(), 40);
        assert_eq!(g1.n_edges(), 40);
        for (_, v, w) in g1.childs(ni(9)) {
            if v == ni(20) {
                assert_abs_diff_eq!(w.trans_prob(), p(1.0));
            } else if v == ni(30) {
                assert!(w.trans_prob().is_zero());
            }
        }
        for (_, v, w) in g1.parents(ni(20)) {
            if v == ni(9) {
                assert_abs_diff_eq!(w.trans_prob(), p(1.0));
            } else if v == ni(19) {
                assert!(w.trans_prob().is_zero());
            }
        }
    }
    #[test]
    fn hmm_crossing_edge_off() {
        let m1 = mock_crossing(false);
        println!("{}", m1);
        let g1 = m1.to_seq_graph().to_phmm(PHMMParams::default());
        println!("{}", g1);
        assert_eq!(g1.n_nodes(), 40);
        assert_eq!(g1.n_edges(), 40);
        for (_, _, w) in g1.childs(ni(9)) {
            assert_abs_diff_eq!(w.trans_prob(), p(0.5));
        }
        for (_, _, w) in g1.parents(ni(20)) {
            assert_abs_diff_eq!(w.trans_prob(), p(0.5));
        }
    }
}
