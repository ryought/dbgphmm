#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ni, sequence_to_string};
    use crate::dbg::mocks::mock_random_with_seq;
    use crate::hmmv2::mocks::{mock_linear_phmm, mock_linear_random_phmm};
    use crate::hmmv2::params::PHMMParams;
    use crate::hmmv2::result::PHMMResultLike;
    use crate::random_seq::generate;
    #[test]
    fn hmmv2_sparse_and_dense() {
        let mut param = PHMMParams::default();
        // param.n_warmup = 2;
        // param.n_active_nodes = 2;
        let (dbg, seq) = mock_random_with_seq(8, 100);
        let phmm = dbg.to_phmm(param);
        println!("{:?}", sequence_to_string(&seq));

        // let r = generate(100, 1);
        let r = &seq[20..60];
        let r2 = &seq[30..70];
        let r1 = phmm.backward(&r);
        let r2 = phmm.backward_sparse(&r2);
        for i in 0..r1.n_emissions() {
            let d = r1.table(i).diff(&r2.table(i));
            println!("{} {}", i, d);
        }
    }
    #[test]
    fn hmmv2_sample_freq() {
        let phmm = mock_linear_phmm(PHMMParams::default());
        // sample
        let h = phmm.sample(10, 4);
        println!("{}", h);
        let r = h.to_sequence();
        let nf_true = h.to_node_freqs(&phmm);
        println!("{:?} {}", nf_true, nf_true.sum());

        // infer
        let o = phmm.run(&r);
        let nf_infer = o.to_node_freqs();
        println!("{:?} {}", nf_infer, nf_infer.sum());

        // infer sparse
        let o = phmm.run_sparse(&r);
        let nf_infer_sparse = o.to_node_freqs();
        println!("{:?} {}", nf_infer_sparse, nf_infer_sparse.sum());
    }
    #[test]
    fn hmmv2_linear_long_node_freq_similarity() {
        // assert that node_freq of
        // * true: determined from the sampling history
        // * dense: estimated with dense calculation
        // * sparse: estimated with sparse calculation
        let phmm = mock_linear_random_phmm(100, 0, PHMMParams::default());
        // sample
        let h = phmm.sample(50, 4);
        let r = h.to_sequence();
        let nf_true = h.to_node_freqs(&phmm);
        println!("{}", nf_true.sum());

        // dense
        let o_dense = phmm.run(&r);
        let nf_dense = o_dense.to_node_freqs();
        println!("{}", nf_dense.sum());

        // sparse
        let o_sparse = phmm.run_sparse(&r);
        let nf_sparse = o_sparse.to_node_freqs();
        println!("{}", nf_sparse.sum());

        // assert node_freq between dense and sparse is similar
        assert_abs_diff_eq!(nf_dense.sum(), nf_sparse.sum(), epsilon = 0.0000001);

        // forward and backward is consistent in dense
        assert_abs_diff_eq!(
            o_dense.to_full_prob_forward(),
            o_dense.to_full_prob_backward(),
            epsilon = 0.1
        );

        // forward and backward is consistent in sparse
        assert_abs_diff_eq!(
            o_sparse.to_full_prob_forward(),
            o_sparse.to_full_prob_backward(),
            epsilon = 0.1
        );

        // full prob of dense and sparse is similar
        assert_abs_diff_eq!(
            o_dense.to_full_prob_forward(),
            o_sparse.to_full_prob_forward(),
            epsilon = 0.1
        );

        // assert all emit probs are similar between dense and sparse
        for (e1, e2) in o_dense.iter_emit_probs().zip(o_sparse.iter_emit_probs()) {
            let diff = e1.diff(&e2);
            println!("{}", diff);
            assert!(diff < 0.000000001);
        }

        for i in 0..r.len() {
            let df = o_dense.forward.table(i);
            let db = o_dense.backward.table(i);
            let sf = o_sparse.forward.table(i);
            let sb = o_sparse.backward.table(i);
            let diff1 = df.diff(&sf);
            let diff2 = db.diff(&sb);
            assert!(diff1 < 0.000000001);
            assert!(diff2 < 0.000000001);
        }
    }
}
