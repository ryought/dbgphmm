//!
//! Hint information
//!
use super::tablev2::{PHMMOutput, MAX_ACTIVE_NODES};
use arrayvec::ArrayVec;
use petgraph::graph::NodeIndex;

/// Hint for a emission sequence
///
/// Define candidate nodes for each emission
///
#[derive(Clone, Debug, PartialEq)]
pub struct Hint(Vec<ArrayVec<NodeIndex, MAX_ACTIVE_NODES>>);

impl Hint {
    /// Constructor from vec of arrayvecs
    ///
    pub fn new(v: Vec<ArrayVec<NodeIndex, MAX_ACTIVE_NODES>>) -> Self {
        Hint(v)
    }
    /// Constructor from vec of vecs
    ///
    pub fn from(v: Vec<Vec<NodeIndex>>) -> Self {
        Hint(
            v.into_iter()
                .map(|nodes| {
                    let mut vec = ArrayVec::new();
                    vec.try_extend_from_slice(&nodes).unwrap();
                    vec
                })
                .collect(),
        )
    }
    /// Get candidate nodes of `emissions[index]`
    ///
    pub fn nodes(&self, index: usize) -> &[NodeIndex] {
        &self.0[index]
    }
    /// Length of the read (emission sequence)
    ///
    pub fn len(&self) -> usize {
        self.0.len()
    }
}

impl PHMMOutput {
    ///
    /// Create Hint of this emission sequence
    ///
    pub fn to_hint(&self, n_active_nodes: usize) -> Hint {
        let ret = self
            .iter_emit_probs()
            .skip(1)
            .map(|state_probs| state_probs.top_nodes(n_active_nodes))
            .collect();
        Hint::new(ret)
    }
}
