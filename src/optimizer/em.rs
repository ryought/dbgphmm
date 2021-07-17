//! EM algorithm
use crate::compressed_dbg::CompressedDBG;
use crate::hmm::base::{PHMMLayer, PHMM};
use crate::hmm::fdbg::FCDbgPHMM;
use crate::hmm::params::PHMMParams;
use rayon::prelude::*;

/// EM optimization with
/// start from initial state
/// -> run forward algorithm to
pub fn optimize_freq_by_em(
    cdbg: &CompressedDBG,
    reads: &[Vec<u8>],
    param: PHMMParams,
    init_freqs: &[f64],
    n_iter: u64,
) {
    let mut freqs = init_freqs.to_vec();

    for i in 0..n_iter {
        println!("freqs={:?}", freqs);

        let phmm = FCDbgPHMM::new(cdbg, freqs);
        let layers: Vec<PHMMLayer> = reads
            .par_iter()
            .map(|read| {
                let f = phmm.forward(&param, read);
                let b = phmm.backward(&param, read);
                let state_prob = phmm.state_prob(&f, &b);
                let ret: PHMMLayer = state_prob.into_iter().sum();
                ret
            })
            .collect();
        let layer_sum: PHMMLayer = layers.into_iter().sum();
        freqs = layer_sum.to_freqs();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_em() {
        let seqs = vec![b"ATTCGATCGATTT".to_vec()];
        let (cdbg, copy_nums) = CompressedDBG::from_seqs(&seqs, 8);
        let freqs = cdbg.copy_nums_to_freqs(&copy_nums);
        let param = PHMMParams::default();
        optimize_freq_by_em(&cdbg, &seqs, param, &freqs, 10);
    }
}
