//!
//! Calculate Node/Edge (i.e. hidden state/transition) usage frequencies
//! from the result of Forward/Backward.
//!
//! - **Emit prob** (for each state and each emission)
//!     The probability of emitting the single base from the hidden state
//!
//!     EmitProb[i, M/I/D, v] = Pr(x[i] is emitted from state M/I/D of node v)
//!
//! - **State probs** (for each state)
//!     The expected value (the sum of probabilities) of the usage frequency
//!     of each hidden state, that is the sum of emit prob, while emitting
//!     the set of emissions.
//!
//!     StateProb[M/I/D, v] = Pr(the expected usage of state M/I/D of node v)
//!
//! - **Node freq** (for each node)
//!     The expected value of the usage frequency of each node.
//!     There are three hidden states for a single node in PHMM, and a node freq
//!     is the sum of three state freqs for `Match/Ins/Del`.
//!
//!     NodeFreq[v] = Pr(the expected usage of node v)
//!
use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::hint::Hint;
use super::tablev2::{NodeVec, PHMMOutput, PHMMTable, PHMMTables, MAX_ACTIVE_NODES};
use super::trans_table::{EdgeFreqs, InitTransProbs, TransProb, TransProbs};
use crate::common::collection::Bases;
use crate::common::{Freq, ReadCollection, Reads, Seq, Sequence};
use crate::graph::active_nodes::ActiveNodes;
use crate::prob::Prob;
use crate::utils::check_memory_usage;
use petgraph::graph::{EdgeIndex, NodeIndex};
use rayon::prelude::*;
use sparsevec::SparseVec;

///
/// methods to generate PHMMOutput from PHMMModel
///
impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Run forward and backward for the emissions and returns PHMMOutput.
    ///
    pub fn run<X: AsRef<Bases>>(&self, emissions: X) -> PHMMOutput {
        let forward = self.forward(&emissions);
        let backward = self.backward(&emissions);
        PHMMOutput::new(forward, backward)
    }
    ///
    /// Run forward and backward with sparse calculation
    /// for the emissions and returns PHMMOutput.
    ///
    pub fn run_sparse<X: AsRef<Bases>>(&self, emissions: X) -> PHMMOutput {
        let forward = self.forward_sparse(&emissions);
        let backward = self.backward_sparse(&emissions);
        PHMMOutput::new(forward, backward)
    }
}

///
/// methods for calculating probabilities
/// with PHMMModel and Reads
///
impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
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
    /// calculate the full probability `P(R)` using rayon parallel calculation.
    ///
    pub fn to_full_prob_parallel<T>(&self, seqs: T) -> Prob
    where
        T: IntoParallelIterator,
        T::Item: Seq,
    {
        seqs.into_par_iter()
            .map(|seq| {
                let read = seq.as_ref();
                let forward = self.forward(read);
                let backward = self.backward(read);
                let o = PHMMOutput::new(forward, backward);
                o.to_full_prob_forward()
            })
            .product()
    }
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
    ///
    /// calculate the full probability `P(R)` using rayon parallel calculation and hint information
    /// (active nodes)
    ///
    /// This function does not run backward. Running forward is enough to calculate the full
    /// probability P(R|G).
    ///
    pub fn to_full_prob_par_with_hint<S>(&self, seqs_and_hints: &[(S, Hint)]) -> Prob
    where
        S: Seq,
    {
        seqs_and_hints
            .into_par_iter()
            .map(|(seq, hint)| {
                let forward = self.forward_with_hint(seq.as_ref(), hint);
                forward.full_prob()
            })
            .product()
    }
    ///
    /// Append hint information in parallel
    ///
    pub fn to_hints_parallel<T>(&self, seqs: T) -> Vec<Hint>
    where
        T: IntoParallelIterator,
        T::Item: Seq,
    {
        seqs.into_par_iter()
            .map(|seq| {
                let hint = self
                    .run_sparse(seq.as_ref())
                    .to_hint(self.param.n_active_nodes);
                hint
            })
            .collect()
    }
}

//
// For hidden states / nodes
//

/// Frequency (f64) assigned to each nodes
pub type NodeFreqs = SparseVec<Freq, NodeIndex, MAX_ACTIVE_NODES>;

/// utility function to convert Vec<Prob> into Vec<f64> by calling Prob::to_log_value
fn from_prob_to_f64_vec(ps: Vec<Prob>) -> Vec<f64> {
    ps.into_iter().map(|p| p.to_log_value()).collect()
}

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
    pub fn iter_emit_probs<'a>(&'a self) -> impl Iterator<Item = PHMMTable> + 'a {
        let n = self.forward.n_emissions();
        let p = self.to_full_prob_forward();
        (0..=n).map(move |i| self.to_emit_probs(i))
    }
    /// Calculate the expected value of the usage frequency of each hidden states
    /// by summing the emit probs of each states for all emissions.
    ///
    pub fn to_state_probs(&self) -> PHMMTable {
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
    ///
    /// Create hint (ActiveNodes list for each bases)
    ///
    pub fn to_hint(&self, n_active_nodes: usize) -> Hint {
        let ret = self
            .iter_emit_probs()
            .skip(1)
            .map(|state_probs| {
                ActiveNodes::Only(state_probs.top_nodes(n_active_nodes).as_slice().to_owned())
            })
            .collect();
        Hint::new(ret)
    }
}

//
// For transitions / edges
//

impl PHMMOutput {
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
    pub fn to_edge_and_init_freqs<N: PHMMNode, E: PHMMEdge>(
        &self,
        phmm: &PHMMModel<N, E>,
        emissions: &[u8],
    ) -> (EdgeFreqs, NodeFreqs) {
        assert_eq!(emissions.len(), self.forward.n_emissions());
        assert_eq!(emissions.len(), self.backward.n_emissions());

        let mut edge_freqs = EdgeFreqs::new(phmm.n_edges(), 0.0, true);
        let mut init_freqs = NodeFreqs::new(phmm.n_nodes(), 0.0, true);

        for i in 0..=self.forward.n_emissions() {
            let (tp, ip) = self.to_trans_and_init_probs(phmm, emissions, i);
            for (e, _, _, _) in phmm.edges() {
                edge_freqs[e] += tp[e].sum().to_value();
            }
            for (v, _) in phmm.nodes() {
                init_freqs[v] += ip[v].sum().to_value();
            }
        }

        (edge_freqs, init_freqs)
    }
    ///
    /// wrapper of to_edge_and_init_freqs
    ///
    pub fn to_edge_freqs<N: PHMMNode, E: PHMMEdge>(
        &self,
        phmm: &PHMMModel<N, E>,
        emissions: &[u8],
    ) -> EdgeFreqs {
        let (ef, _) = self.to_edge_and_init_freqs(phmm, emissions);
        ef
    }
    ///
    /// wrapper of PHMMResult.to_trans_and_init_probs
    ///
    pub fn to_trans_probs<N: PHMMNode, E: PHMMEdge>(
        &self,
        phmm: &PHMMModel<N, E>,
        emissions: &[u8],
        i: usize,
    ) -> TransProbs {
        let (tp, _) = self.to_trans_and_init_probs(phmm, emissions, i);
        tp
    }
    /// Calculate the expected value of the usage frequency of each edges
    ///
    /// ```text
    /// P(pi[i]=k and pi[i+1]=l)
    ///  = P(emits x[:i] and pi=k) p_trans(k,l) p_emit(l, x[i]) P(emits x[i+1:] | pi=l)
    ///  = f[i][k] * p_trans * p_emit * b[i+1][l]
    /// ```
    ///
    /// range of i: 0 <= i <= n
    ///
    pub fn to_trans_and_init_probs<N: PHMMNode, E: PHMMEdge>(
        &self,
        phmm: &PHMMModel<N, E>,
        emissions: &[u8],
        i: usize,
    ) -> (TransProbs, InitTransProbs) {
        let n = emissions.len();
        assert_eq!(n, self.forward.n_emissions());
        assert_eq!(n, self.backward.n_emissions());
        assert!(i <= n);

        let mut t: TransProbs = TransProbs::new(phmm.n_edges(), TransProb::zero(), true);
        let mut t0: InitTransProbs = InitTransProbs::new(phmm.n_nodes(), TransProb::zero(), true);

        let param = &phmm.param;
        let p = self.to_full_prob_forward();

        // get table references
        let fi0 = self.forward.table_merged(i);
        let bi2 = self.backward.table_merged(i + 1);
        let bi1 = self.backward.table_merged(i);

        // (1) fill TransProbs
        for (e, k, l, ew) in phmm.edges() {
            let p_trans = ew.trans_prob();

            // to m (normal state)
            if i < n {
                let p_emit = phmm.p_match_emit(l, emissions[i]);
                t[e].mm = fi0.m[k] * p_trans * param.p_MM * p_emit * bi2.m[l] / p;
                t[e].im = fi0.i[k] * p_trans * param.p_IM * p_emit * bi2.m[l] / p;
                t[e].dm = fi0.d[k] * p_trans * param.p_DM * p_emit * bi2.m[l] / p;
            }

            // to d (silent state)
            t[e].md = fi0.m[k] * p_trans * param.p_MD * bi1.d[l] / p;
            t[e].id = fi0.i[k] * p_trans * param.p_ID * bi1.d[l] / p;
            t[e].dd = fi0.d[k] * p_trans * param.p_DD * bi1.d[l] / p;
        }

        // (2) fill InitTransProbs
        for (v, vw) in phmm.nodes() {
            let p_init = vw.init_prob();

            // to m (normal state)
            if i < n {
                let p_emit = phmm.p_match_emit(v, emissions[i]);
                t0[v].mm = fi0.mb * p_init * param.p_MM * p_emit * bi2.m[v] / p;
                t0[v].im = fi0.ib * p_init * param.p_IM * p_emit * bi2.m[v] / p;
            }

            // to d (silent state)
            t0[v].md = fi0.mb * p_init * param.p_MD * bi1.d[v] / p;
            t0[v].id = fi0.ib * p_init * param.p_ID * bi1.d[v] / p;
        }

        (t, t0)
    }
}

//
// For inspection of emissions
//

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// PHMMModel version of `Dbg::show_mapping_summary`
    ///
    fn inspect<F: Fn(NodeIndex) -> String>(
        &self,
        output: &PHMMOutput,
        emissions: &[u8],
        node_info: F,
    ) {
        assert_eq!(output.n_emissions(), emissions.len());
        // for (i, &emission) in emissions.into_iter().enumerate() {}
        for (i, state_prob) in output.iter_emit_probs().skip(1).enumerate() {
            // println!("i={} {}", i, state_prob);
            println!(
                "{}\t{}\t{}",
                i,
                emissions[i] as char,
                state_prob.to_summary_string(&node_info)
            );
        }
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
        let eps: Vec<PHMMTable> = o.iter_emit_probs().collect();
        for ep in eps.iter() {
            println!("{}", ep);
        }
        assert_abs_diff_eq!(eps[5].m[ni(7)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(eps[4].m[ni(6)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(eps[3].m[ni(5)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(eps[2].m[ni(4)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(eps[1].m[ni(3)], p(1.0), epsilon = 0.00001);
        assert_abs_diff_eq!(eps[0].mb, p(1.0), epsilon = 0.00001);

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
        for i in 0..=es.len() {
            let (tp, ip) = o.to_trans_and_init_probs(&phmm, es, i);
            println!("{}", i);
            phmm.draw_edge_vec(&tp);
            phmm.draw_node_vec(&ip);
        }
        assert!(o
            .to_trans_probs(&phmm, es, 0)
            .iter()
            .all(|(_, t)| t.sum().is_zero()));
        assert_abs_diff_eq!(
            o.to_trans_probs(&phmm, es, 1)[ei(3)].mm,
            p(1.0),
            epsilon = 0.00001
        );
        assert_abs_diff_eq!(
            o.to_trans_probs(&phmm, es, 2)[ei(4)].mm,
            p(1.0),
            epsilon = 0.00001
        );
        assert_abs_diff_eq!(
            o.to_trans_probs(&phmm, es, 3)[ei(5)].mm,
            p(1.0),
            epsilon = 0.00001
        );
        assert_abs_diff_eq!(
            o.to_trans_probs(&phmm, es, 4)[ei(6)].mm,
            p(1.0),
            epsilon = 0.00001
        );
        assert!(o
            .to_trans_probs(&phmm, es, 5)
            .iter()
            .all(|(_, t)| t.sum().is_zero()));

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
        let tps: Vec<TransProbs> = (1..=es.len())
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
    fn hmm_node_freq_inspect() {
        let phmm = mock_linear_phmm(PHMMParams::default());
        let read = b"ATTCGTCGT";
        let output = phmm.run(read);
        phmm.inspect(&output, read, |_| format!(""));
    }
}
