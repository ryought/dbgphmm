//!
//! HMMv2 test
//!
//!
pub mod dbg;

#[cfg(test)]
mod tests {
    use crate::common::{ni, sequence_to_string};
    use crate::e2e;
    use crate::genome;
    use crate::hmmv2::mocks::{mock_linear_phmm, mock_linear_random_phmm};
    use crate::hmmv2::params::PHMMParams;
    use crate::prelude::*;
    use crate::random_seq::generate;
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
}
