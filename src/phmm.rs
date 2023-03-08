//!
//! Profile hidden Markov model (PHMM) for sequence graph
//!
//! * Infer
//!     Estimating output probability or emitting states given emissions
//! * Sample
//!     Sampling emissions (reads) from the model
//!
//! # Usage
//!
//! Define PHMM that is petgraph::DiGraph in which each node corresponds to a base
//!
pub mod infer;
pub mod params;
pub mod sample;

// re-export
pub use params::PHMMParams;

use crate::common::{CopyNum, NULL_BASE};
use crate::prob::Prob;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

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
pub struct PHMM {
    param: PHMMParams,
    graph: DiGraph<PHMMNode, PHMMEdge>,
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
