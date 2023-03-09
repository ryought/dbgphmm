//!
//! Profile hidden Markov model (PHMM) for sequence graph
//!
//! * Infer
//!     Estimating output probability or emitting states given emissions
//!     i.e. (PHMM, emissions) -> information about emissions
//! * Sample
//!     Sampling emissions (reads) from the model
//!     i.e. PHMM -> History (sequence of emission and originating state)
//!
//! # Usage
//!
//! Define PHMM that is petgraph::DiGraph in which each node corresponds to a base
//!
pub mod infer;
pub mod nodesubset;
pub use nodesubset::NodeSubset;
pub mod params;
pub mod sample;

// re-export
pub use params::PHMMParams;

use crate::common::NULL_BASE;
use crate::prob::Prob;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph_algos::iterators::{ChildEdges, EdgesIterator, NodesIterator, ParentEdges};

///
/// Node in PHMM
///
/// # Attributes
///
/// * `emission`
///     A base for emission of this node
/// * `init_prob`
///     Initial probability from Begin node into this node
///
pub struct PHMMNode {
    ///
    /// A base for emission of this node
    ///
    emission: u8,
    ///
    /// Initial probability from Begin node into this node
    ///
    init_prob: Prob,
}

///
/// Edge in PHMM
///
/// # Attributes
///
/// * `trans_prob`
///     Transition probability from source node into target node of this edge
///
pub struct PHMMEdge {
    ///
    /// Transition probability from source node into target node of this edge
    ///
    trans_prob: Prob,
}

///
/// Profile HMM
///
/// * graph: DiGraph<PHMMNode, PHMMEdge>
/// * param: PHMMParams
///
pub struct PHMM {
    param: PHMMParams,
    graph: DiGraph<PHMMNode, PHMMEdge>,
}

///
/// PHMM States
///
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum State {
    Match(NodeIndex),
    Ins(NodeIndex),
    Del(NodeIndex),
    MatchBegin,
    InsBegin,
    End,
}

///
/// PHMM emission
///
/// Del states emit no base, so it will be marked by `Empty`.
///
#[derive(Debug, Copy, Clone)]
pub enum Emission {
    Base(u8),
    Empty,
}

//
//
// implementations
//
//

//
// Node
//

impl PHMMNode {
    ///
    /// Constructor `PHMMNode::new(init_prob, emission)`
    ///
    pub fn new(init_prob: Prob, emission: u8) -> PHMMNode {
        PHMMNode {
            init_prob,
            emission,
        }
    }
    fn emission(&self) -> u8 {
        self.emission
    }
    fn is_emittable(&self) -> bool {
        self.emission() != NULL_BASE
    }
    fn init_prob(&self) -> Prob {
        self.init_prob
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

//
// Edge
//

impl PHMMEdge {
    pub fn new(trans_prob: Prob) -> PHMMEdge {
        PHMMEdge { trans_prob }
    }
    fn trans_prob(&self) -> Prob {
        self.trans_prob
    }
}

impl std::fmt::Display for PHMMEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "p={}", self.trans_prob)
    }
}

//
// Model
//

impl PHMM {
    pub fn param(&self) -> PHMMParams {
        self.param
    }
    pub fn graph(&self) -> &DiGraph<PHMMNode, PHMMEdge> {
        &self.graph
    }
    pub fn n_nodes(&self) -> usize {
        &self.graph.node_count()
    }
    ///
    /// Emission of the node
    ///
    pub fn emission(&self, node: NodeIndex) -> u8 {
        self.graph.node_weight(node).unwrap().emission()
    }
    ///
    /// Init prob of the node
    ///
    pub fn init_prob(&self, node: NodeIndex) -> Prob {
        self.graph.node_weight(node).unwrap().init_prob()
    }
    ///
    /// Trans prob of the edge
    ///
    pub fn trans_prob(&self, edge: EdgeIndex) -> Prob {
        self.graph.edge_weight(edge).unwrap().trans_prob()
    }
    ///
    /// Iterator of all nodes in the graph
    ///
    /// Item is (node: NodeIndex, node_weight: &PHMMNode)
    ///
    pub fn nodes(&self) -> NodesIterator<PHMMNode> {
        NodesIterator::new(&self.graph)
    }
    ///
    /// Iterator of all edges in the graph
    ///
    /// Item is (edge: EdgeIndex, source: NodeIndex, target: NodeIndex, edge_weight: &PHMMEdge)
    ///
    pub fn edges(&self) -> EdgesIterator<PHMMEdge> {
        EdgesIterator::new(&self.graph)
    }
    ///
    /// Iterator of childs of the node
    ///
    /// Item is (edge: EdgeIndex, child: NodeIndex, edge_weight: &PHMMEdge)
    ///
    pub fn childs(&self, node: NodeIndex) -> ChildEdges<PHMMEdge> {
        ChildEdges::new(&self.graph, node)
    }
    ///
    /// Iterator of parents of the node
    ///
    /// Item is (edge: EdgeIndex, parent: NodeIndex, edge_weight: &PHMMEdge)
    ///
    pub fn parents(&self, node: NodeIndex) -> ParentEdges<PHMMEdge> {
        ParentEdges::new(&self.graph, node)
    }
}

impl PHMM {
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

//
// State
//

impl State {
    ///
    /// Convert State (either Match/Ins/Del(NodeIndex) MatchBegin/InsBegin/End)
    /// into the wrapped NodeIndex.
    /// If state is begin/end, it returns None.
    ///
    pub fn to_node_index(&self) -> Option<NodeIndex> {
        match self {
            State::Match(v) | State::Ins(v) | State::Del(v) => Some(*v),
            _ => None,
        }
    }
}

impl std::fmt::Display for State {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            State::Match(node) => write!(f, "M(n{})", node.index()),
            State::Ins(node) => write!(f, "I(n{})", node.index()),
            State::Del(node) => write!(f, "D(n{})", node.index()),
            State::MatchBegin => write!(f, "MB"),
            State::InsBegin => write!(f, "IB"),
            State::End => write!(f, "E"),
        }
    }
}

//
// Emission
//

impl Emission {
    ///
    /// Check if the emission is actual base or not.
    ///
    pub fn is_base(&self) -> bool {
        match self {
            Emission::Base(_) => true,
            Emission::Empty => false,
        }
    }
}

impl std::fmt::Display for Emission {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Emission::Base(base) => write!(f, "{}", *base as char),
            Emission::Empty => write!(f, "-"),
        }
    }
}
