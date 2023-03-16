use super::common::PHMMModel;
use super::table::{NodeVec, PHMMOutput, PHMMTable, PHMMTables, MAX_ACTIVE_NODES};
use crate::common::{Freq, ReadCollection, Reads, Seq, Sequence};
use crate::prob::Prob;
use petgraph::graph::{EdgeIndex, NodeIndex};
use rayon::prelude::*;
use sparsevec::SparseVec;
// use crate::utils::check_memory_usage;

///
/// methods to generate PHMMOutput from PHMMModel
///
impl PHMMModel {
    ///
    /// Run forward and backward for the emissions and returns PHMMOutput.
    ///
    pub fn run(&self, emissions: &[u8]) -> PHMMOutput {
        let forward = self.forward(emissions);
        let backward = self.backward(emissions);
        PHMMOutput::new(forward, backward)
    }
    ///
    /// Run forward and backward with sparse calculation
    /// for the emissions and returns PHMMOutput.
    ///
    pub fn run_sparse(&self, emissions: &[u8]) -> PHMMOutput {
        let forward = self.forward_sparse(emissions);
        let backward = self.backward_sparse(emissions);
        PHMMOutput::new(forward, backward)
    }
}

///
/// methods for calculating probabilities
/// with PHMMModel and Reads
///
impl PHMMModel {
    ///
    /// calculate node freqs of multiple emission sequences (`Reads`).
    ///
    pub fn to_node_freqs<T>(&self, seqs: T) -> NodeFreqs
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        seqs.into_iter()
            .map(|seq| {
                let read = seq.as_ref();
                let forward = self.forward(read);
                let backward = self.backward(read);
                let o = PHMMOutput::new(forward, backward);
                o.to_node_freqs()
            })
            .sum()
    }
    ///
    ///
    ///
    pub fn to_full_prob<T>(&self, seqs: T) -> Prob
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        seqs.into_iter()
            .map(|seq| {
                let read = seq.as_ref();
                let forward = self.forward(read);
                let backward = self.backward(read);
                let o = PHMMOutput::new(forward, backward);
                o.to_full_prob_forward()
            })
            .product()
    }
    /// Calculate full prob using sparse
    ///
    pub fn to_full_prob_sparse<T>(&self, seqs: T) -> Prob
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        seqs.into_iter()
            .map(|seq| {
                let read = seq.as_ref();
                let forward = self.forward_sparse(read);
                forward.full_prob()
            })
            .product()
    }
    /// Calculate full prob using sparse
    ///
    pub fn to_full_prob_sparse_backward<T>(&self, seqs: T) -> Prob
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        seqs.into_iter()
            .map(|seq| {
                let read = seq.as_ref();
                self.backward_sparse(read).full_prob()
            })
            .product()
    }
}

//
// For hidden states / nodes
//

/// The probability of emitting the emission from the hidden state.
pub type EmitProbs = Vec<PHMMTable>;

/// Probability assigned to each hidden states
pub type StateProbs = PHMMTable;

/// Frequency (f64) assigned to each nodes
pub type NodeFreqs = SparseVec<f64, NodeIndex, MAX_ACTIVE_NODES>;

impl PHMMOutput {
    /// Calculate the probability that the hidden states (that is (type, node))
    /// emits the i-th state.
    ///
    /// ## Description
    ///
    /// HMM state probability
    ///
    /// P(pi=s)
    ///  = P(emits x[:i] and pi=s) P(emits x[i:] | pi=s)
    ///  = f[i] + b[i]
    ///
    /// ## Old Description
    ///
    /// `freq[i][t_k]`
    /// = P(base `x[i]` is emitted from node `t_k`)
    /// = `(ft_i[k] * bt_i+1[k]) / P(x)`
    ///
    /// * `t` is a type of state, either Match, Ins, Del
    /// * `k` is a node index
    /// * `P(x)` is the full probability of emissions
    ///
    pub fn iter_emit_probs<'a>(&'a self) -> impl Iterator<Item = StateProbs> + 'a {
        let n = self.forward.n_emissions();
        let p = self.to_full_prob_forward();
        (0..=n).map(move |i| {
            let f = self.forward.table_merged(i);
            let b = self.backward.table_merged(i);
            (f * b) / p
        })
    }
    /// Calculate the expected value of the usage frequency of each hidden states
    /// by summing the emit probs of each states for all emissions.
    ///
    pub fn to_state_probs(&self) -> StateProbs {
        self.iter_emit_probs().sum()
    }
    /// Calculate the expected value of the usage frequency of each nodes
    /// by summing the emit probs of M/I/D states for each node.
    /// `f[v]` = (How many times the hidden state `M_v, I_v, D_v` was used to emit the whole
    /// emissions?)
    ///
    pub fn to_node_freqs(&self) -> NodeFreqs {
        // v is NodeVec<Prob>
        let state_probs = self.to_state_probs();
        let v = state_probs.to_nodevec();
        // f is NodeVec<Freq>
        let mut f: NodeFreqs = NodeFreqs::new(state_probs.n_nodes(), 0.0, true);
        for (node, p) in v.iter() {
            f[node] = p.to_value();
        }
        f
    }
}

//
// Tests
//

#[cfg(test)]
mod tests {}
