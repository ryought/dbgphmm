//!
//! Hint information
//!
use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::table::{PHMMOutput, MAX_ACTIVE_NODES};
use crate::common::{ReadCollection, Seq};
use crate::prob::Prob;
use crate::utils::progress_common_style;
use arrayvec::ArrayVec;
use indicatif::ParallelProgressIterator;
use petgraph::graph::NodeIndex;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

///
/// Mapping information for reads
///
pub type Mappings = Vec<Hint>;

/// Hint for a emission sequence
///
/// Define candidate nodes for each emission
///
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
pub struct Hint {
    nodes: Vec<ArrayVec<NodeIndex, MAX_ACTIVE_NODES>>,
    probs: Vec<ArrayVec<Prob, MAX_ACTIVE_NODES>>,
}

impl Hint {
    /// Constructor from vec of arrayvecs
    ///
    pub fn new(nodes: Vec<ArrayVec<NodeIndex, MAX_ACTIVE_NODES>>) -> Self {
        Hint {
            nodes,
            probs: vec![],
        }
    }
    ///
    /// Constructor from vec of vecs
    ///
    pub fn from(v: Vec<Vec<NodeIndex>>) -> Self {
        Hint {
            nodes: v
                .into_iter()
                .map(|nodes| {
                    let mut vec = ArrayVec::new();
                    // TODO
                    vec.try_extend_from_slice(&nodes).unwrap();
                    vec
                })
                .collect(),
            probs: vec![],
        }
    }
    ///
    pub fn to_inner(self) -> Vec<ArrayVec<NodeIndex, MAX_ACTIVE_NODES>> {
        self.nodes
    }
    ///
    pub fn inner(&self) -> &[ArrayVec<NodeIndex, MAX_ACTIVE_NODES>] {
        &self.nodes
    }
    /// Get candidate nodes of `emissions[index]`
    ///
    pub fn nodes(&self, index: usize) -> &[NodeIndex] {
        &self.nodes[index]
    }
    /// Length of the read (emission sequence)
    ///
    pub fn len(&self) -> usize {
        self.nodes.len()
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

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    ///
    pub fn generate_mappings<S: Seq>(
        &self,
        reads: &ReadCollection<S>,
        mappings: Option<&Mappings>,
    ) -> Mappings {
        reads
            .par_iter()
            .enumerate()
            .progress_with_style(progress_common_style())
            .map(|(i, seq)| {
                let output = if let Some(mappings) = mappings {
                    self.run_with_hint(seq.as_ref(), &mappings[i])
                } else {
                    self.run_sparse(seq.as_ref())
                };
                output.to_hint(self.param.n_active_nodes)
            })
            .collect()
    }
}
