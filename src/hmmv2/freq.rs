//!
//! Calculate Node/Edge (i.e. hidden state/transition) usage frequencies
//! from the result of Forward/Backward.
//!
use super::common::{PHMMEdge, PHMMModel, PHMMNode};

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    pub fn to_node_freq() {}
    pub fn to_edge_freq() {}
}
