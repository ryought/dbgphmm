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
use super::result::{PHMMResult, PHMMResultLike, PHMMResultSparse};
use super::table::PHMMTable;
use super::table_ref::PHMMTableRef;
use super::trans_table::{EdgeFreqs, TransProb, TransProbs};
use crate::common::{Freq, Reads, Sequence};
use crate::prob::Prob;
use crate::vector::{DenseStorage, EdgeVec, NodeVec, Storage};
use petgraph::graph::{EdgeIndex, NodeIndex};

/// Struct for storing `PHMMResultLike` for forward and backward.
///
#[derive(Debug, Clone)]
pub struct PHMMOutput<R: PHMMResultLike> {
    /// PHMMResult for forward run
    pub forward: R,
    /// PHMMResult for backward run
    pub backward: R,
}

/// Constructors
impl<R: PHMMResultLike> PHMMOutput<R> {
    fn new(forward: R, backward: R) -> Self {
        // TODO check forward is created by phmm.forward()
        PHMMOutput { forward, backward }
    }
}

///
/// methods to generate PHMMOutput from PHMMModel
///
impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Run forward and backward for the emissions and returns PHMMOutput.
    ///
    pub fn run(&self, emissions: &[u8]) -> PHMMOutput<PHMMResult> {
        let forward = self.forward(emissions);
        let backward = self.backward(emissions);
        PHMMOutput::new(forward, backward)
    }
    ///
    /// Run forward and backward with sparse calculation
    /// for the emissions and returns PHMMOutput.
    ///
    pub fn run_sparse(&self, emissions: &[u8]) -> PHMMOutput<PHMMResultSparse> {
        let forward = self.forward_sparse(emissions);
        let backward = self.backward_sparse(emissions);
        PHMMOutput::new(forward, backward)
    }
}

///
/// methods for calculating probabilities
/// with PHMMModel and Reads
///
impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    /// calculate node freqs of multiple emission sequences (`Reads`).
    pub fn to_node_freqs(&self, r: &Reads) -> NodeFreqs {
        r.reads
            .iter()
            .map(|read| {
                let forward = self.forward(read);
                let backward = self.backward(read);
                let o = PHMMOutput::new(forward, backward);
                o.to_node_freqs()
            })
            .sum()
    }
    /// calculate edge freqs of multiple emission sequences (`Reads`).
    pub fn to_edge_freqs(&self, r: &Reads) -> EdgeFreqs {
        r.reads
            .iter()
            .map(|read| {
                let forward = self.forward(read);
                let backward = self.backward(read);
                let o = PHMMOutput::new(forward, backward);
                o.to_edge_freqs(self, read)
            })
            .sum()
    }
}

//
// For full probability
//

impl<R: PHMMResultLike> PHMMOutput<R> {
    /// Calculate the full probability `P(x)` of the given emission `x`
    /// from **forward** result.
    ///
    /// ```text
    /// P(x) = fe_n-1 = P(emits x[0],...,x[n-1] and now in `e` (end state))
    /// ```
    ///
    pub fn to_full_prob_forward(&self) -> Prob {
        match self.forward.last_table() {
            PHMMTableRef::Dense(t) => t.e,
            PHMMTableRef::Sparse(t) => t.e,
        }
    }
    /// Calculate the full probability `P(x)` of the given emission `x`
    /// from **backward** result.
    ///
    /// ```text
    /// P(x) = bm_0[b] = P(emits x[0:] | starts from m_b)
    /// ```
    ///
    pub fn to_full_prob_backward(&self) -> Prob {
        match self.backward.first_table() {
            PHMMTableRef::Dense(t) => t.mb,
            PHMMTableRef::Sparse(t) => t.mb,
        }
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

impl<R: PHMMResultLike> PHMMOutput<R> {
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
    pub fn iter_emit_probs<'a>(&'a self) -> impl Iterator<Item = StateProbs> + 'a {
        let n = self.forward.n_emissions();
        let p = self.to_full_prob_forward();
        (0..n).map(move |i| {
            let f = self.forward.table(i);
            let b = if i + 1 < n {
                self.backward.table(i + 1)
            } else {
                self.backward.init_table()
            };
            (&f * &b) / p
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
        let mut f: NodeFreqs = NodeFreqs::new(state_probs.n_nodes(), 0.0);
        for (node, p) in v.iter() {
            f[node] = p.to_value();
        }
        f
    }
}

//
// For transitions / edges
//

impl<R: PHMMResultLike> PHMMOutput<R> {
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
    pub fn to_edge_freqs<N: PHMMNode, E: PHMMEdge>(
        &self,
        phmm: &PHMMModel<N, E>,
        emissions: &[u8],
    ) -> EdgeFreqs {
        assert_eq!(emissions.len(), self.forward.n_emissions());
        assert_eq!(emissions.len(), self.backward.n_emissions());

        let mut freq: EdgeFreqs = EdgeFreqs::new(phmm.n_edges(), 0.0);

        for i in 0..self.forward.n_emissions() {
            let tp = self.to_trans_probs(phmm, emissions, i);
            for (e, _, _, _) in phmm.edges() {
                freq[e] += tp[e].sum().to_value();
            }
        }

        freq
    }
    /// Calculate the expected value of the usage frequency of each edges
    /// TBW
    pub fn to_trans_probs<N: PHMMNode, E: PHMMEdge>(
        &self,
        phmm: &PHMMModel<N, E>,
        emissions: &[u8],
        i: usize,
    ) -> TransProbs {
        assert_eq!(emissions.len(), self.forward.n_emissions());
        assert_eq!(emissions.len(), self.backward.n_emissions());

        let mut t: TransProbs = TransProbs::new(phmm.n_edges(), TransProb::zero());

        let param = &phmm.param;
        let fi0 = self.forward.table(i);
        let p = self.to_full_prob_forward();
        let n = self.forward.n_emissions();

        // to m (normal state)
        let bi2 = if i + 2 < n {
            self.backward.table(i + 2)
        } else {
            self.backward.init_table()
        };
        let bi1 = if i + 1 < n {
            self.backward.table(i + 1)
        } else {
            self.backward.init_table()
        };

        if i + 1 < n {
            for (e, k, l, ew) in phmm.edges() {
                let p_emit = phmm.p_match_emit(l, emissions[i + 1]);
                let p_trans = ew.trans_prob();
                t[e].mm = fi0.m(k) * p_trans * param.p_MM * p_emit * bi2.m(l) / p;
                t[e].im = fi0.i(k) * p_trans * param.p_IM * p_emit * bi2.m(l) / p;
                t[e].dm = fi0.d(k) * p_trans * param.p_DM * p_emit * bi2.m(l) / p;
            }
        }

        // to d (silent state)
        if i < n {
            for (e, k, l, weight) in phmm.edges() {
                let p_trans = weight.trans_prob();
                t[e].md = fi0.m(k) * p_trans * param.p_MD * bi1.d(l) / p;
                t[e].id = fi0.i(k) * p_trans * param.p_ID * bi1.d(l) / p;
                t[e].dd = fi0.d(k) * p_trans * param.p_DD * bi1.d(l) / p;
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
    use crate::hmmv2::mocks::mock_linear_phmm;
    use crate::hmmv2::params::PHMMParams;
    use crate::prob::p;
    use crate::vector::DenseStorage;
    #[test]
    fn hmm_freq_mock_linear_zero_error_full_prob() {
        let phmm = mock_linear_phmm(PHMMParams::zero_error());
        let o = phmm.run(b"CGATC");
        assert_abs_diff_eq!(
            o.to_full_prob_forward(),
            o.to_full_prob_backward(),
            epsilon = 0.0000001,
        );
    }
    #[test]
    fn hmm_freq_mock_linear_zero_error_node_freqs() {
        let phmm = mock_linear_phmm(PHMMParams::zero_error());
        let o = phmm.run(b"CGATC");
        let eps: Vec<StateProbs> = o.iter_emit_probs().collect();
        for ep in eps.iter() {
            println!("{}", ep);
        }
        assert_abs_diff_eq!(eps[4].m[ni(7)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(eps[3].m[ni(6)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(eps[2].m[ni(5)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(eps[1].m[ni(4)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(eps[0].m[ni(3)], p(1.0), epsilon = 0.00001);

        let sps = o.to_state_probs();
        println!("{}", sps);
        assert_abs_diff_eq!(sps.m[ni(7)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(sps.m[ni(6)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(sps.m[ni(5)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(sps.m[ni(4)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(sps.m[ni(3)], p(1.0), epsilon = 0.00001);

        let nf = o.to_node_freqs();
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
        let o = phmm.run(b"CGATC");
        let sps = o.to_state_probs();
        let nf = o.to_node_freqs();
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
        let o = phmm.run(es);
        let sps = o.to_state_probs();
        let nf = o.to_node_freqs();
        phmm.draw_node_vec(&nf);
        assert_abs_diff_eq!(
            o.to_full_prob_forward(),
            o.to_full_prob_backward(),
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
        let o = phmm.run(es);

        // (1) trans_probs
        for i in 0..5 {
            let tps = o.to_trans_probs(&phmm, es, i);
            println!("{}", i);
            phmm.draw_edge_vec(&tps);
        }
        assert_abs_diff_eq!(
            o.to_trans_probs(&phmm, es, 0)[ei(3)].mm,
            p(1.0),
            epsilon = 0.00001
        );
        assert_abs_diff_eq!(
            o.to_trans_probs(&phmm, es, 1)[ei(4)].mm,
            p(1.0),
            epsilon = 0.00001
        );
        assert_abs_diff_eq!(
            o.to_trans_probs(&phmm, es, 2)[ei(5)].mm,
            p(1.0),
            epsilon = 0.00001
        );
        assert_abs_diff_eq!(
            o.to_trans_probs(&phmm, es, 3)[ei(6)].mm,
            p(1.0),
            epsilon = 0.00001
        );

        // edge_freqs
        let efs = o.to_edge_freqs(&phmm, es);
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
        let o = phmm.run(es);
        for (i, table) in o.forward.tables.iter().enumerate() {
            println!("rf[{}]\n{}", i, table);
        }
        for (i, table) in o.backward.tables.iter().enumerate() {
            println!("rb[{}]\n{}", i, table);
        }

        // (1) trans_probs
        let tps: Vec<TransProbs> = (0..es.len())
            .map(|i| o.to_trans_probs(&phmm, es, i))
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
        let efs = o.to_edge_freqs(&phmm, es);
        phmm.draw_edge_vec(&efs);
        // all edge will be used so all should have f~1
        for (e, _, _, _) in phmm.edges() {
            assert_abs_diff_eq!(efs[e], 0.99, epsilon = 0.01);
        }
    }
    #[test]
    fn hmm_to_node_freqs() {
        let phmm = mock_linear_phmm(PHMMParams::default());
    }
}
