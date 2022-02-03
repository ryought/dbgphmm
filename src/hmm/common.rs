//!
//! Definition of Node-centric PHMM
//!
use super::params::PHMMParams;
use crate::common::CopyNum;
use crate::prob::Prob;
use petgraph::dot::Dot;
use petgraph::graph::DiGraph;
use petgraph::graph::Edges;
use petgraph::graph::NodeReferences;
pub use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::{EdgeRef, IntoNodeReferences};
use petgraph::Directed;
pub use petgraph::Direction;

//
// traits
//
pub trait PHMMNode {
    fn emission(&self) -> u8;
    fn is_emittable(&self) -> bool {
        self.emission() != b'N'
    }
    // fn copy_num(&self) -> CopyNum;
    fn init_prob(&self) -> Prob;
}

pub trait PHMMEdge {
    // fn copy_num(&self) -> CopyNum;
    fn trans_prob(&self) -> Prob;
}

pub struct PHMMModel<N: PHMMNode, E: PHMMEdge> {
    pub param: PHMMParams,
    pub graph: DiGraph<N, E>,
}

pub type PGraph = DiGraph<PNode, PEdge>;
pub type PModel = PHMMModel<PNode, PEdge>;

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// create iterator of all nodes
    /// Item of the iterator is `(NodeIndex, &N)`
    ///
    pub fn nodes(&self) -> Nodes<N> {
        Nodes {
            nodes: (&self.graph).node_references(),
        }
    }
    ///
    /// create iterator of all child edges of the node
    /// Item of the iterator is EdgeReference `er`
    ///
    pub fn childs(&self, v: NodeIndex) -> ChildEdges<E> {
        ChildEdges {
            edges: self.graph.edges_directed(v, Direction::Outgoing),
        }
    }
    ///
    /// create iterator of all parent edges of the node
    /// Item of the iterator is EdgeReference `er`
    ///
    pub fn parents(&self, v: NodeIndex) -> ParentEdges<E> {
        ParentEdges {
            edges: self.graph.edges_directed(v, Direction::Incoming),
        }
    }
}

pub struct Nodes<'a, N: 'a> {
    nodes: NodeReferences<'a, N>,
}

impl<'a, N> Iterator for Nodes<'a, N> {
    type Item = (NodeIndex, &'a N);
    fn next(&mut self) -> Option<Self::Item> {
        self.nodes.next()
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
