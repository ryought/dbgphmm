//! EM algorithm
use crate::compressed_dbg::CompressedDBG;
use crate::hmm::base::{PHMMLayer, PHMM};
use crate::hmm::cdbg::CDbgPHMM;
use crate::hmm::fdbg::FCDbgPHMM;
use crate::hmm::params::PHMMParams;
use crate::optimizer::annealer::Annealer;
use crate::optimizer::freq::FreqState;
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
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

///
pub fn optimize_copy_nums_by_em(
    cdbg: &CompressedDBG,
    reads: &[Vec<u8>],
    param: PHMMParams,
    init_copy_nums: &[u32],
) {
    let mut copy_nums = init_copy_nums.to_vec();

    // E-step: copy_nums -> freqs
    let phmm = CDbgPHMM::new(cdbg, copy_nums);
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
    let freqs: Vec<f64> = layer_sum.to_freqs().iter().map(|f| f / 20.0).collect();
    println!("E: freqs={:?}", freqs);

    // M-step: freqs -> copy_nums
    let s = FreqState::init(&cdbg, &freqs);
    let a = Annealer::new(1.0, 0.8);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
    a.run_with_log(&mut rng, s, 100);
}

///
pub fn optimize_copy_nums_by_em_with_true(
    cdbg: &CompressedDBG,
    reads: &[Vec<u8>],
    param: PHMMParams,
    copy_nums_true: &[u32],
    copy_nums_read: &[u32],
) {
    let freqs: Vec<f64> = copy_nums_read.iter().map(|&cn| cn as f64 / 20.0).collect();

    println!("{:?}", copy_nums_true);
    println!("{:?}", copy_nums_read);
    println!("{:?}", freqs);

    // let s = FreqState::new(&cdbg, &freqs, copy_nums_true);
    let s = FreqState::init(&cdbg, &freqs);
    let a = Annealer::new(1.0, 0.8);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
    a.run_with_log(&mut rng, s, 100);
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
