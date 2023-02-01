#[cfg(test)]
mod tests {
    use crate::common::{ni, sequence_to_string};
    use crate::dbg::mocks::mock_random_with_seq;
    use crate::e2e;
    use crate::genome;
    use crate::hmmv2::mocks::{mock_linear_phmm, mock_linear_random_phmm};
    use crate::hmmv2::params::PHMMParams;
    use crate::hmmv2::result::PHMMResultLike;
    use crate::prelude::*;
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
    fn hmmv2_incorrect_sparse_calculation() {
        let dataset = e2e::generate_tandem_repeat_1kbp();
        // dataset.show_reads_with_genome();

        let k = 40;
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k, dataset.genome());
        let phmm = dbg.to_phmm(PHMMParams::uniform(0.001));
        let read = &dataset.reads()[0];
        let output = phmm.run_sparse(read.seq());
        println!("p_f={}", output.to_full_prob_forward());
        println!("p_b={}", output.to_full_prob_backward());
    }
}
