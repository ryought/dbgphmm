//!
//! Calculate Node/Edge (i.e. hidden state/transition) usage frequencies
//! from the result of Forward/Backward.
//!
//! - **Emit prob** (for each state and each emission)
//!     The probability of emitting the single base from the hidden state
//!
//! - **State probs** (for each state)
//!     The expected value (the sum of probabilities) of the usage frequency
//!     of each hidden state, that is the sum of emit prob, while emitting
//!     the set of emissions.
//!
//! - **Node freq** (for each node)
//!     The expected value of the usage frequency of each node.
//!     There are three hidden states for a single node in PHMM, and a node freq
//!     is the sum of three state freqs for `Match/Ins/Del`.
//!
use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::table::{PHMMResult, PHMMTable};
use crate::common::Freq;
use crate::prob::Prob;
use crate::vector::{DenseStorage, EdgeVec, NodeVec, Storage};
use petgraph::graph::NodeIndex;

/// The probability of emitting the emission from the hidden state.
pub type EmitProbs = Vec<PHMMTable<DenseStorage<Prob>>>;

/// Probability assigned to each hidden states
pub type StateProbs = PHMMTable<DenseStorage<Prob>>;

/// Frequency (f64) assigned to each nodes
pub type NodeFreqs = NodeVec<DenseStorage<Freq>>;

/// Frequency (f64) assigned to each edges
pub type EdgeFreqs = EdgeVec<DenseStorage<Freq>>;

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    /// Calculate the full probability `P(x)` of the given emission `x`
    /// from **forward** result.
    ///
    /// ```text
    /// P(x) = fe_n-1 = P(emits x[0],...,x[n-1] and now in `e` (end state))
    /// ```
    ///
    pub fn to_full_prob_forward<S>(&self, forward: &PHMMResult<S>) -> Prob
    where
        S: Storage<Item = Prob>,
    {
        let f = forward.tables.last().unwrap();
        f.e
    }
    /// Calculate the full probability `P(x)` of the given emission `x`
    /// from **backward** result.
    ///
    /// ```text
    /// P(x) = bm_0[b] = P(emits x[0:] | starts from m_b)
    /// ```
    ///
    pub fn to_full_prob_backward<S>(&self, backward: &PHMMResult<S>) -> Prob
    where
        S: Storage<Item = Prob>,
    {
        let b = backward.tables.first().unwrap();
        b.mb
    }
    /// Calculate the probability that the hidden states (that is (type, node))
    /// emits the i-th emission.
    ///
    /// `freq[i][t_k]`
    /// = P(base `x[i]` is emitted from node `t_k`)
    /// = `ft_i[k] * bt_i+1[k]`
    ///
    /// * `t` is a type of state, either Match, Ins, Del
    /// * `k` is a node index
    ///
    pub fn to_emit_probs<S>(&self, forward: &PHMMResult<S>, backward: &PHMMResult<S>) -> EmitProbs
    where
        S: Storage<Item = Prob>,
    {
        let n = forward.n_emissions();
        (0..n)
            .map(|i| {
                let f = &forward.tables[i];
                let b = if i + 1 < n {
                    &backward.tables[i + 1]
                } else {
                    &backward.init_table
                };
                f * b
            })
            .map(|v| v.to_dense())
            .collect()
    }
    /// Calculate the expected value of the usage frequency of each hidden states
    /// by summing the emit probs of each states for all emissions.
    ///
    pub fn to_state_probs<S>(&self, forward: &PHMMResult<S>, backward: &PHMMResult<S>) -> StateProbs
    where
        S: Storage<Item = Prob>,
    {
        self.to_emit_probs(forward, backward).into_iter().sum()
    }
    /// Calculate the expected value of the usage frequency of each nodes
    /// by summing the emit probs of M/I/D states for each node.
    /// `f[v]` = (How many times the hidden state `M_v, I_v, D_v` was used to emit the whole
    /// emissions?)
    ///
    pub fn to_node_freqs<S>(&self, state_probs: &StateProbs) -> NodeFreqs
    where
        S: Storage<Item = Prob>,
    {
        let n = state_probs.n_nodes();
        let mut f: NodeFreqs = NodeFreqs::new(n, 0.0);
        for i in 0..n {
            let v = NodeIndex::new(i);
            f[v] = (state_probs.m[v] + state_probs.i[v] + state_probs.d[v]).to_value();
        }
        f
    }
    /// Calculate the expected value of the usage frequency of each edges
    ///
    pub fn to_edge_freqs<S>(&self, forward: &PHMMResult<S>, backward: &PHMMResult<S>) -> EdgeFreqs
    where
        S: Storage<Item = Prob>,
    {
        unimplemented!();
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
    use crate::vector::DenseStorage;
    #[test]
    fn hmm_freq_mock_linear() {
        let phmm = mock_linear()
            .to_seq_graph()
            .to_phmm(PHMMParams::high_error());
        let rf: PHMMResult<DenseStorage<Prob>> = phmm.forward(b"CGATC");
        let rb: PHMMResult<DenseStorage<Prob>> = phmm.backward(b"CGATC");
        // phmm.to_node_freq(&rf, &rb);
    }
}
