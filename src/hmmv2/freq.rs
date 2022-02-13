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
use super::result::PHMMResult;
use super::table::PHMMTable;
use super::trans_table::{EdgeFreqs, TransProb, TransProbs};
use crate::common::Freq;
use crate::prob::Prob;
use crate::vector::{DenseStorage, EdgeVec, NodeVec, Storage};
use petgraph::graph::{EdgeIndex, NodeIndex};

//
// For full probability
//

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
}

//
// For hidden states / nodes
//

/// The probability of emitting the emission from the hidden state.
pub type EmitProbs = Vec<PHMMTable<DenseStorage<Prob>>>;

/// Probability assigned to each hidden states
pub type StateProbs = PHMMTable<DenseStorage<Prob>>;

/// Frequency (f64) assigned to each nodes
pub type NodeFreqs = NodeVec<DenseStorage<Freq>>;

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    /// Calculate the probability that the hidden states (that is (type, node))
    /// emits the i-th emission.
    ///
    /// `freq[i][t_k]`
    /// = P(base `x[i]` is emitted from node `t_k`)
    /// = `(ft_i[k] * bt_i+1[k]) / P(x)`
    ///
    /// * `t` is a type of state, either Match, Ins, Del
    /// * `k` is a node index
    /// * `P(x)` is the full probability of emissions
    ///
    pub fn to_emit_probs<S>(&self, forward: &PHMMResult<S>, backward: &PHMMResult<S>) -> EmitProbs
    where
        S: Storage<Item = Prob>,
    {
        let n = forward.n_emissions();
        let p = self.to_full_prob_forward(forward);
        (0..n)
            .map(|i| {
                let f = &forward.tables[i];
                let b = if i + 1 < n {
                    &backward.tables[i + 1]
                } else {
                    &backward.init_table
                };
                (f * b) / p
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
        // TODO to_emit_probs can be an iterator (storeing all temp vector is unnecessary).
        self.to_emit_probs(forward, backward).into_iter().sum()
    }
    /// Calculate the expected value of the usage frequency of each nodes
    /// by summing the emit probs of M/I/D states for each node.
    /// `f[v]` = (How many times the hidden state `M_v, I_v, D_v` was used to emit the whole
    /// emissions?)
    ///
    /// TODO this does not depend on self, so should move to StateProbs method?
    ///
    pub fn to_node_freqs(&self, state_probs: &StateProbs) -> NodeFreqs {
        let n = state_probs.n_nodes();
        let mut f: NodeFreqs = NodeFreqs::new(n, 0.0);
        for i in 0..n {
            let v = NodeIndex::new(i);
            f[v] = (state_probs.m[v] + state_probs.i[v] + state_probs.d[v]).to_value();
        }
        f
    }
}

//
// For transitions / edges
//

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    /// Calculate the expected value of the usage frequency of each edges
    ///
    /// `freq[i][e]`
    /// = P(use an edge `e` after emitting x[i])
    ///
    /// For a Normal state `N` and silent state `S`,
    ///
    /// * `e = (k: N or S, l: N)` transition
    ///     `freq[i][e] = (f_i[k] * a_kl * e_l(x[i+1]) * b_i+2[l]) / P(x)`
    ///
    /// * `e = (k: N or S, l: S)` transition
    ///     `freq[i][e] = (f_i[k] * a_kl * b_i+1[l]) / P(x)`
    ///
    pub fn to_edge_freqs<S>(
        &self,
        forward: &PHMMResult<S>,
        backward: &PHMMResult<S>,
        emissions: &[u8],
    ) -> EdgeFreqs
    where
        S: Storage<Item = Prob>,
    {
        assert_eq!(emissions.len(), forward.n_emissions());
        assert_eq!(emissions.len(), backward.n_emissions());

        let mut freq: EdgeFreqs = EdgeFreqs::new(self.n_edges(), 0.0);

        for i in 0..forward.n_emissions() {
            let tp = self.to_trans_probs(forward, backward, emissions, i);
            for (e, _, _, _) in self.edges() {
                freq[e] += tp[e].sum().to_value();
            }
        }

        freq
    }
    /// Calculate the expected value of the usage frequency of each edges
    /// TBW
    pub fn to_trans_probs<S>(
        &self,
        forward: &PHMMResult<S>,
        backward: &PHMMResult<S>,
        emissions: &[u8],
        i: usize,
    ) -> TransProbs
    where
        S: Storage<Item = Prob>,
    {
        assert_eq!(emissions.len(), forward.n_emissions());
        assert_eq!(emissions.len(), backward.n_emissions());

        let mut t: TransProbs = TransProbs::new(self.n_edges(), TransProb::zero());

        let param = &self.param;
        let fi0 = &forward.tables[i];
        let p = self.to_full_prob_forward(forward);

        // to m (normal state)
        let bi2 = if i + 2 < forward.n_emissions() {
            &backward.tables[i + 2]
        } else {
            &backward.init_table
        };
        let bi1 = if i + 1 < forward.n_emissions() {
            &backward.tables[i + 1]
        } else {
            &backward.init_table
        };

        if i + 1 < forward.n_emissions() {
            for (e, k, l, ew) in self.edges() {
                let p_emit = self.p_match_emit(l, emissions[i + 1]);
                let p_trans = ew.trans_prob();
                t[e].mm = fi0.m[k] * p_trans * param.p_MM * p_emit * bi2.m[l] / p;
                t[e].im = fi0.i[k] * p_trans * param.p_IM * p_emit * bi2.m[l] / p;
                t[e].dm = fi0.d[k] * p_trans * param.p_DM * p_emit * bi2.m[l] / p;
            }
        }

        // to d (silent state)
        if i < forward.n_emissions() {
            for (e, k, l, weight) in self.edges() {
                let p_trans = weight.trans_prob();
                t[e].md = fi0.m[k] * p_trans * param.p_MD * bi1.d[l] / p;
                t[e].id = fi0.i[k] * p_trans * param.p_ID * bi1.d[l] / p;
                t[e].dd = fi0.d[k] * p_trans * param.p_DD * bi1.d[l] / p;
            }
        }
        t
    }
}

//
// Tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};
    use crate::hmm::params::PHMMParams;
    use crate::hmmv2::mocks::mock_linear_phmm;
    use crate::prob::p;
    use crate::vector::DenseStorage;
    #[test]
    fn hmm_freq_mock_linear_zero_error_full_prob() {
        let phmm = mock_linear_phmm(PHMMParams::zero_error());
        let rf: PHMMResult<DenseStorage<Prob>> = phmm.forward(b"CGATC");
        let rb: PHMMResult<DenseStorage<Prob>> = phmm.backward(b"CGATC");
        assert_abs_diff_eq!(
            phmm.to_full_prob_forward(&rf),
            phmm.to_full_prob_backward(&rb),
            epsilon = 0.0000001,
        );
    }
    #[test]
    fn hmm_freq_mock_linear_zero_error_node_freqs() {
        let phmm = mock_linear_phmm(PHMMParams::zero_error());
        let rf: PHMMResult<DenseStorage<Prob>> = phmm.forward(b"CGATC");
        let rb: PHMMResult<DenseStorage<Prob>> = phmm.backward(b"CGATC");
        let eps = phmm.to_emit_probs(&rf, &rb);
        for ep in eps.iter() {
            println!("{}", ep);
        }
        assert_abs_diff_eq!(eps[4].m[ni(7)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(eps[3].m[ni(6)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(eps[2].m[ni(5)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(eps[1].m[ni(4)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(eps[0].m[ni(3)], p(1.0), epsilon = 0.00001);

        let sps = phmm.to_state_probs(&rf, &rb);
        println!("{}", sps);
        assert_abs_diff_eq!(sps.m[ni(7)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(sps.m[ni(6)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(sps.m[ni(5)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(sps.m[ni(4)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(sps.m[ni(3)], p(1.0), epsilon = 0.00001);

        let nf = phmm.to_node_freqs(&sps);
        println!("{:?}", nf);
        assert_abs_diff_eq!(nf[ni(7)], 1.0, epsilon = 0.00001);
        assert_abs_diff_eq!(nf[ni(6)], 1.0, epsilon = 0.00001);
        assert_abs_diff_eq!(nf[ni(5)], 1.0, epsilon = 0.00001);
        assert_abs_diff_eq!(nf[ni(4)], 1.0, epsilon = 0.00001);
        assert_abs_diff_eq!(nf[ni(3)], 1.0, epsilon = 0.00001);
        assert_abs_diff_eq!(nf[ni(9)], 0.0, epsilon = 0.00001);
        assert_abs_diff_eq!(nf[ni(8)], 0.0, epsilon = 0.00001);
        assert_abs_diff_eq!(nf[ni(2)], 0.0, epsilon = 0.00001);
        assert_abs_diff_eq!(nf[ni(1)], 0.0, epsilon = 0.00001);
        assert_abs_diff_eq!(nf[ni(0)], 0.0, epsilon = 0.00001);
    }
    #[test]
    fn hmm_freq_mock_linear_high_error_node_freqs() {
        let phmm = mock_linear_phmm(PHMMParams::default());
        let rf: PHMMResult<DenseStorage<Prob>> = phmm.forward(b"CGATC");
        let rb: PHMMResult<DenseStorage<Prob>> = phmm.backward(b"CGATC");
        let sps = phmm.to_state_probs(&rf, &rb);
        let nf = phmm.to_node_freqs(&sps);
        phmm.draw_node_vec(&nf);
        assert!(nf[ni(0)] < 0.01); // -
        assert!(nf[ni(1)] < 0.01); // -
        assert!(nf[ni(2)] < 0.01); // -
        assert!(nf[ni(3)] > 0.98); // C
        assert!(nf[ni(4)] > 0.98); // G
        assert!(nf[ni(5)] > 0.98); // A
        assert!(nf[ni(6)] > 0.98); // T
        assert!(nf[ni(7)] > 0.98); // C
        assert!(nf[ni(8)] < 0.01); // -
        assert!(nf[ni(9)] < 0.01); // -
    }
    #[test]
    fn hmm_freq_mock_linear_error_node_freqs_del() {
        let phmm = mock_linear_phmm(PHMMParams::default());
        // orig: b"ATTCGATCGT";
        let es = b"ATTCGTCGT"; // have 1 deletion
        let rf: PHMMResult<DenseStorage<Prob>> = phmm.forward(es);
        let rb: PHMMResult<DenseStorage<Prob>> = phmm.backward(es);
        let sps = phmm.to_state_probs(&rf, &rb);
        let nf = phmm.to_node_freqs(&sps);
        phmm.draw_node_vec(&nf);
        assert_abs_diff_eq!(
            phmm.to_full_prob_forward(&rf),
            phmm.to_full_prob_backward(&rb),
            epsilon = 0.00001,
        );
        for (v, _) in phmm.nodes() {
            println!("{}", nf[v]);
            assert_abs_diff_eq!(nf[v], 1.0, epsilon = 0.01);
        }
    }
    #[test]
    fn hmm_freq_mock_linear_zero_error_trans_probs() {
        let phmm = mock_linear_phmm(PHMMParams::zero_error());
        let es = b"CGATC";
        let rf: PHMMResult<DenseStorage<Prob>> = phmm.forward(es);
        let rb: PHMMResult<DenseStorage<Prob>> = phmm.backward(es);

        // (1) trans_probs
        for i in 0..5 {
            let tps = phmm.to_trans_probs(&rf, &rb, es, i);
            println!("{}", i);
            phmm.draw_edge_vec(&tps);
        }
        assert_abs_diff_eq!(
            phmm.to_trans_probs(&rf, &rb, es, 0)[ei(3)].mm,
            p(1.0),
            epsilon = 0.00001
        );
        assert_abs_diff_eq!(
            phmm.to_trans_probs(&rf, &rb, es, 1)[ei(4)].mm,
            p(1.0),
            epsilon = 0.00001
        );
        assert_abs_diff_eq!(
            phmm.to_trans_probs(&rf, &rb, es, 2)[ei(5)].mm,
            p(1.0),
            epsilon = 0.00001
        );
        assert_abs_diff_eq!(
            phmm.to_trans_probs(&rf, &rb, es, 3)[ei(6)].mm,
            p(1.0),
            epsilon = 0.00001
        );

        // edge_freqs
        let efs = phmm.to_edge_freqs(&rf, &rb, es);
        phmm.draw_edge_vec(&efs);
        assert!(efs[ei(0)] < 0.0001);
        assert!(efs[ei(1)] < 0.0001);
        assert!(efs[ei(2)] < 0.0001);
        assert!(efs[ei(3)] > 0.9999); // C -> G
        assert!(efs[ei(4)] > 0.9999); // G -> A
        assert!(efs[ei(5)] > 0.9999); // A -> T
        assert!(efs[ei(6)] > 0.9999); // T -> C
        assert!(efs[ei(7)] < 0.0001);
        assert!(efs[ei(8)] < 0.0001);
    }
    #[test]
    fn hmm_freq_mock_linear_high_error_trans_probs() {
        let phmm = mock_linear_phmm(PHMMParams::default());
        // let es = b"ATTCGATCGT";
        let es = b"ATTCGTCGT";
        let rf: PHMMResult<DenseStorage<Prob>> = phmm.forward(es);
        let rb: PHMMResult<DenseStorage<Prob>> = phmm.backward(es);
        for (i, table) in rf.tables.iter().enumerate() {
            println!("rf[{}]\n{}", i, table);
        }
        for (i, table) in rb.tables.iter().enumerate() {
            println!("rb[{}]\n{}", i, table);
        }

        // (1) trans_probs
        let tps: Vec<TransProbs> = (0..es.len())
            .map(|i| phmm.to_trans_probs(&rf, &rb, es, i))
            .collect();
        for (i, t) in tps.iter().enumerate() {
            println!("{}", i);
            phmm.draw_edge_vec(t);
        }
        assert!(tps[0][ei(0)].mm.to_value() > 0.9);
        assert!(tps[1][ei(1)].mm.to_value() > 0.9);
        assert!(tps[2][ei(2)].mm.to_value() > 0.9);
        assert!(tps[3][ei(3)].mm.to_value() > 0.9);
        assert!(tps[4][ei(4)].md.to_value() > 0.9);
        assert!(tps[4][ei(5)].dm.to_value() > 0.9);
        assert!(tps[5][ei(6)].mm.to_value() > 0.9);
        assert!(tps[6][ei(7)].mm.to_value() > 0.9);
        assert!(tps[7][ei(8)].mm.to_value() > 0.9);

        // (2) edge_freqs
        println!("edge_freqs");
        let efs = phmm.to_edge_freqs(&rf, &rb, es);
        phmm.draw_edge_vec(&efs);
        // all edge will be used so all should have f~1
        for (e, _, _, _) in phmm.edges() {
            assert_abs_diff_eq!(efs[e], 0.99, epsilon = 0.01);
        }
    }
}
