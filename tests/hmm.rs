//!
//! test of hmm
//!
#[macro_use]
extern crate approx;

use dbgphmm::common::{ni, sequence_to_string};
use dbgphmm::dbg::mocks::mock_random_with_seq;
use dbgphmm::hmmv2::common::PModel;
use dbgphmm::hmmv2::mocks::{mock_linear_phmm, mock_linear_random_phmm};
use dbgphmm::hmmv2::params::PHMMParams;
use dbgphmm::hmmv2::result::PHMMResultLike;
use dbgphmm::random_seq::generate;

fn check_node_freq_similarity(phmm: &PModel) {
    // assert that node_freq of
    // * true: determined from the sampling history
    // * dense: estimated with dense calculation
    // * sparse: estimated with sparse calculation
    // sample
    let h = phmm.sample(500, 4);
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

    // forward/backward table is similar regardless of dense/sparse calculation
    for i in 0..r.len() {
        let df = o_dense.forward.table(i);
        let db = o_dense.backward.table(i);
        let sf = o_sparse.forward.table(i);
        let sb = o_sparse.backward.table(i);
        let diff1 = df.diff(&sf);
        let diff2 = db.diff(&sb);
        println!("forward_dense-forward_sparse = {}", diff1);
        println!("backward_dense-backward_sparse = {}", diff1);
        assert!(diff1 < 0.000000001);
        assert!(diff2 < 0.000000001);
    }

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

    // assert node_freq between dense and sparse is similar
    assert_abs_diff_eq!(nf_dense.sum(), nf_sparse.sum(), epsilon = 0.0000001);
}

#[test]
#[ignore]
fn hmmv2_linear_long_node_freq_similarity() {
    let phmm = mock_linear_random_phmm(1000, 2, PHMMParams::default());
    check_node_freq_similarity(&phmm);
}
