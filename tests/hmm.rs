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
use dbgphmm::hmmv2::sample::History;
use dbgphmm::random_seq::generate;
use itertools::izip;

fn check_node_freq_similarity(phmm: &PModel, h: &History) {
    // assert that node_freq of
    // * true: determined from the sampling history
    // * dense: estimated with dense calculation
    // * sparse: estimated with sparse calculation
    // sample
    println!("sampling reads..");
    let r = h.to_sequence();
    println!("read created");
    let nf_true = h.to_node_freqs(&phmm);
    println!("{}", nf_true.sum());

    // dense
    println!("running dense...");
    let o_dense = phmm.run(&r);
    let nf_dense = o_dense.to_node_freqs();
    println!("{}", nf_dense.sum());

    // sparse
    println!("running sparse...");
    let o_sparse = phmm.run_sparse(&r);
    let nf_sparse = o_sparse.to_node_freqs();
    println!("{}", nf_sparse.sum());

    /*
    for (v, _) in phmm.nodes() {
        println!("{:?}: dense - sparse={}", v, nf_dense[v] - nf_sparse[v]);
    }
    */

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
        // println!("forward_dense-forward_sparse = {}", diff1);
        // println!("backward_dense-backward_sparse = {}", diff1);
        assert!(diff1 < 0.000000001);
        assert!(diff2 < 0.000000001);
    }

    // full prob of dense and sparse is similar
    assert_abs_diff_eq!(
        o_dense.to_full_prob_forward(),
        o_sparse.to_full_prob_forward(),
        epsilon = 0.1
    );
    println!("p_f[dense, forward]={}", o_dense.to_full_prob_forward());
    println!("p_f[sparse, forward]={}", o_sparse.to_full_prob_forward());
    println!("p_f[dense, backward]={}", o_dense.to_full_prob_backward());
    println!("p_f[sparse, backward]={}", o_sparse.to_full_prob_backward());

    // assert all emit probs are similar between dense and sparse
    let epd: Vec<_> = o_dense.iter_emit_probs().collect();
    let eps: Vec<_> = o_sparse.iter_emit_probs().collect();
    for (i, (e1, e2)) in izip!(epd, eps).enumerate() {
        let diff = e1.diff(&e2);
        /*
        if diff > 0.000000001 {
            println!("e1={}", e1);
            println!("e2={}", e2);
            println!("m");
            e1.m.show_diff(&e2.m);
            println!("i");
            e1.i.show_diff(&e2.i);
            println!("d");
            e1.d.show_diff(&e2.d);

            let df = o_dense.forward.table(i);
            let db = o_dense.backward.table(i + 1);
            let sf = o_sparse.forward.table(i);
            let sb = o_sparse.backward.table(i + 1);
            println!("df(is_dense={})", df.is_dense());
            println!("db(is_dense={})", db.is_dense());
            println!("sf(is_dense={})", sf.is_dense());
            println!("sb(is_dense={})", sb.is_dense());
            println!("diff of f {}", df.diff(&sf));
            println!("diff of f {}", sf.diff(&df));
            println!("diff of b {}", db.diff(&sb));
            println!("diff of b {}", sb.diff(&db));
        }
        */
        println!("diff={}", diff);
        assert!(diff < 0.001);
    }

    // assert node_freq between dense and sparse is similar
    assert_abs_diff_eq!(nf_dense.sum(), nf_sparse.sum(), epsilon = 0.01);
}

#[test]
#[ignore]
fn hmmv2_linear_long_node_freq_similarity() {
    let mut param = PHMMParams::default();
    param.n_warmup = 40;
    // param.n_warmup = 20;
    // param.n_active_nodes = 16;
    // param.n_active_nodes = 32;
    let phmm = mock_linear_random_phmm(1000, 2, param);
    println!("phmm created");
    let h = phmm.sample(500, 4);
    // let h = phmm.sample(500, 2);
    // let h = phmm.sample(500, 3);
    // println!("{}", h);
    println!("sample created");
    check_node_freq_similarity(&phmm, &h);
}
