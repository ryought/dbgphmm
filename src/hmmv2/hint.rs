//!
//! Hint information
//!
use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::table::{PHMMOutput, MAX_ACTIVE_NODES};
use crate::common::{ReadCollection, Seq};
use arrayvec::ArrayVec;
use indicatif::{ParallelProgressIterator, ProgressStyle};
use petgraph::graph::NodeIndex;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Hint for a emission sequence
///
/// Define candidate nodes for each emission
///
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
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
    ///
    pub fn to_inner(self) -> Vec<ArrayVec<NodeIndex, MAX_ACTIVE_NODES>> {
        self.0
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

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Append hint information
    ///
    /// # Todos
    ///
    /// * make parallel by `into_par_iter`
    ///
    pub fn append_hints<S: Seq>(
        &self,
        reads: ReadCollection<S>,
        parallel: bool,
    ) -> ReadCollection<S> {
        let style = ProgressStyle::with_template(
            "[{elapsed_precise}/{eta_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-");

        let hints = if parallel {
            reads
                .par_iter()
                .progress_with_style(style)
                .map(|seq| {
                    self.run_sparse(seq.as_ref())
                        .to_hint(self.param.n_active_nodes)
                })
                .collect()
        } else {
            reads
                .into_iter()
                .enumerate()
                .map(|(i, seq)| {
                    println!("sparse... #{}", i);
                    self.run_sparse(seq.as_ref())
                        .to_hint(self.param.n_active_nodes)
                })
                .collect()
        };
        ReadCollection::from_with_hint(reads.reads, hints)
    }
}
