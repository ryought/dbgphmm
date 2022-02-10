//!
//! Definition of Node-centric PHMM
//!
use crate::common::CopyNum;
use crate::hmm::params::PHMMParams;
use crate::prob::Prob;
use crate::vector::{EdgeVec, NodeVec, Storage};
use petgraph::dot::Dot;
use petgraph::graph::DiGraph;
pub use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::graph::{EdgeReferences, Edges, NodeReferences};
use petgraph::visit::{EdgeRef, IntoNodeReferences};
use petgraph::Directed;
pub use petgraph::Direction;

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
pub trait PHMMNode {
    ///
    /// Emission base assigned to this HMM node.
    fn emission(&self) -> u8;
    ///
    /// This node has a valid emission or not.
    fn is_emittable(&self) -> bool {
        self.emission() != b'N'
    }
    ///
    /// Initial transition probability from Begin state to the node.
    fn init_prob(&self) -> Prob;
    // TODO
    // fn copy_num(&self) -> CopyNum;
}

///
/// attribute of PHMM Edge
///
/// * `trans_prob(&self) -> Prob`
///     Transition probability from `source` to `target`.
///
pub trait PHMMEdge {
    ///
    /// Transition probability from the source node to the target node
    /// of this edge.
    fn trans_prob(&self) -> Prob;
    // TODO
    // fn copy_num(&self) -> CopyNum;
}

/// Profile HMM model
///
pub struct PHMMModel<N: PHMMNode, E: PHMMEdge> {
    pub param: PHMMParams,
    pub graph: DiGraph<N, E>,
}

pub type PGraph = DiGraph<PNode, PEdge>;
pub type PModel = PHMMModel<PNode, PEdge>;

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// create iterator of all nodes
    /// Item of the iterator is `(NodeIndex, &N)`.
    ///
    pub fn nodes(&self) -> NodesIterator<N> {
        NodesIterator {
            nodes: (&self.graph).node_references(),
        }
    }
    ///
    /// create iterator of all nodes
    /// Item of the iterator is
    /// `(EdgeIndex, NodeIndex of source, NodeIndex of target, &E)`.
    ///
    pub fn edges(&self) -> EdgesIterator<E> {
        EdgesIterator {
            edges: (&self.graph).edge_references(),
        }
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
    pub fn childs(&self, node: NodeIndex) -> ChildEdges<E> {
        ChildEdges {
            edges: self.graph.edges_directed(node, Direction::Outgoing),
        }
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
    pub fn parents(&self, node: NodeIndex) -> ParentEdges<E> {
        ParentEdges {
            edges: self.graph.edges_directed(node, Direction::Incoming),
        }
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
    /// Get a node emission
    ///
    pub fn emission(&self, v: NodeIndex) -> u8 {
        self.graph.node_weight(v).unwrap().emission()
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

/// Iterator struct for `nodes()`
/// (with `NodeReferences` of petgraph)
///
pub struct NodesIterator<'a, N: 'a> {
    nodes: NodeReferences<'a, N>,
}

impl<'a, N> Iterator for NodesIterator<'a, N> {
    type Item = (NodeIndex, &'a N);
    fn next(&mut self) -> Option<Self::Item> {
        self.nodes.next()
    }
}

/// Iterator struct for `edges()`
/// (with `NodeReferences` of petgraph)
///
pub struct EdgesIterator<'a, E: 'a> {
    edges: EdgeReferences<'a, E>,
}

impl<'a, E> Iterator for EdgesIterator<'a, E> {
    type Item = (EdgeIndex, NodeIndex, NodeIndex, &'a E);
    fn next(&mut self) -> Option<Self::Item> {
        // extract edge reference
        match self.edges.next() {
            Some(er) => Some((er.id(), er.source(), er.target(), er.weight())),
            None => None,
        }
    }
}

pub struct ChildEdges<'a, E: 'a> {
    edges: Edges<'a, E, Directed>,
}

impl<'a, E> Iterator for ChildEdges<'a, E> {
    type Item = (EdgeIndex, NodeIndex, &'a E);
    fn next(&mut self) -> Option<Self::Item> {
        // edge reference
        match self.edges.next() {
            // er.source() = the given node
            // er.target() = child
            Some(er) => Some((er.id(), er.target(), er.weight())),
            None => None,
        }
    }
}

pub struct ParentEdges<'a, E: 'a> {
    edges: Edges<'a, E, Directed>,
}

impl<'a, E> Iterator for ParentEdges<'a, E> {
    type Item = (EdgeIndex, NodeIndex, &'a E);
    fn next(&mut self) -> Option<Self::Item> {
        // edge reference
        match self.edges.next() {
            // er.source() = parent
            // er.target() = the given node
            Some(er) => Some((er.id(), er.source(), er.weight())),
            None => None,
        }
    }
}

impl<N, E> std::fmt::Display for PHMMModel<N, E>
where
    N: PHMMNode + std::fmt::Display,
    E: PHMMEdge + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.graph, &[]))
    }
}

//
//
// PNode
//
//

#[derive(Debug, Copy, Clone)]
pub struct PNode {
    ///
    copy_num: CopyNum,
    /// initial transition probability
    ///  = (copy_num) / (sum of copy_nums of all nodes)
    init_prob: Prob,
    ///
    is_emittable: bool,
    ///
    emission: u8,
}

impl PNode {
    pub fn new(copy_num: CopyNum, init_prob: Prob, is_emittable: bool, emission: u8) -> PNode {
        PNode {
            copy_num,
            init_prob,
            is_emittable,
            emission,
        }
    }
}

impl PHMMNode for PNode {
    fn emission(&self) -> u8 {
        self.emission
    }
    fn is_emittable(&self) -> bool {
        self.is_emittable
    }
    fn init_prob(&self) -> Prob {
        self.init_prob
    }
}

impl std::fmt::Display for PNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.is_emittable {
            write!(f, "{} (p={})", self.emission as char, self.init_prob)
        } else {
            write!(f, "not_emittable (p={})", self.init_prob)
        }
    }
}

//
//
// PEdge
//
//

#[derive(Debug, Copy, Clone)]
pub struct PEdge {
    ///
    trans_prob: Prob,
    ///
    copy_num: Option<CopyNum>,
}

impl PEdge {
    pub fn new(trans_prob: Prob) -> PEdge {
        PEdge {
            trans_prob,
            copy_num: None,
        }
    }
}

impl PHMMEdge for PEdge {
    fn trans_prob(&self) -> Prob {
        self.trans_prob
    }
}

impl std::fmt::Display for PEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "p={}", self.trans_prob)
    }
}

//
// Visualizers
//
impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    /// println wrapper of NodeVec, with attributes of each nodes
    pub fn draw_node_vec<S>(&self, nv: &NodeVec<S>)
    where
        S: Storage,
        S::Item: std::fmt::Display,
    {
        for (node, _w) in self.nodes() {
            println!("{:?}({})\t{}", node, self.emission(node) as char, nv[node]);
        }
    }
    /// println wrapper of EdgeVec, with attributes of each edges
    pub fn draw_edge_vec<S>(&self, ev: &EdgeVec<S>)
    where
        S: Storage,
        S::Item: std::fmt::Display,
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
    use crate::graph::mocks::mock_linear;
    use crate::hmm::params::PHMMParams;

    #[test]
    fn phmmmodels_basic_ops() {
        let phmm = mock_linear()
            .to_seq_graph()
            .to_phmm(PHMMParams::high_error());
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
}
