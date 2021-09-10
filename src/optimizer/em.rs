//! EM algorithm
use crate::compressed_dbg::CompressedDBG;
use crate::hmm::base::{PHMMLayer, PHMM};
use crate::hmm::cdbg::CDbgPHMM;
use crate::hmm::fdbg::FCDbgPHMM;
use crate::hmm::params::PHMMParams;
use crate::optimizer::annealer::Annealer;
use crate::optimizer::bestfreq::BestFreqState;
use crate::optimizer::freq::FreqState;
use crate::optimizer::grad::GradientDescent;
use crate::prob::Prob;
use log::{info, warn};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

/*
 * schedulers
 */
pub trait DepthScheduler {
    fn depth(&self, iteration: u64) -> f64;
}

#[derive(Debug, Clone)]
pub struct ConstantDepth {
    depth: f64,
}

impl ConstantDepth {
    pub fn new(depth: f64) -> ConstantDepth {
        ConstantDepth { depth }
    }
}

impl DepthScheduler for ConstantDepth {
    fn depth(&self, _: u64) -> f64 {
        self.depth
    }
}

#[derive(Debug, Clone)]
pub struct LinearGradientDepth {
    init_depth: f64,
    final_depth: f64,
    n_iter: u64,
}

impl LinearGradientDepth {
    pub fn new(init_depth: f64, final_depth: f64, n_iter: u64) -> LinearGradientDepth {
        LinearGradientDepth {
            init_depth,
            final_depth,
            n_iter,
        }
    }
}

impl DepthScheduler for LinearGradientDepth {
    fn depth(&self, iteration: u64) -> f64 {
        let d0 = self.init_depth as f64;
        let d = self.final_depth as f64;
        let t = iteration as f64 / (self.n_iter as f64 - 1.0);
        d0 + ((d - d0) * t)
    }
}

/// EM iterative optimization on freq f64 space.
/// experimental function
pub fn optimize_freq_by_em(
    cdbg: &CompressedDBG,
    reads: &[Vec<u8>],
    param: PHMMParams,
    init_freqs: &[f64],
    n_iter: u64,
) {
    let mut freqs = init_freqs.to_vec();

    for i in 0..n_iter {
        let phmm = FCDbgPHMM::new(cdbg, freqs.clone());
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
        let p_total: Prob = reads
            .par_iter()
            .map(|read| phmm.forward_prob(&param, read))
            .product();

        println!(
            "{}\t{:.16}\t{:.16}\t{}\tnull\t{:?}",
            i,
            1.0,
            p_total.to_log_value(),
            "",
            freqs
        );

        freqs = layer_sum.to_freqs();
    }
}

/// EM iterative optimization on copy_nums space
/// In E-step, freqs of kmers is computed with PHMM with given copy_nums
/// Next in M-step, copy_nums is updated to that best fits with freqs
pub fn optimize_copy_nums_by_em<T: DepthScheduler>(
    cdbg: &CompressedDBG,
    reads: &[Vec<u8>],
    param: PHMMParams,
    init_copy_nums: &[u32],
    depth_scheduler: &T,
    n_iter: u64,
) {
    let mut copy_nums = init_copy_nums.to_vec();

    for i in 0..n_iter {
        // E-step: copy_nums -> freqs
        let phmm = CDbgPHMM::new(cdbg, copy_nums.clone());
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
        let depth = depth_scheduler.depth(i);
        let freqs: Vec<f64> = layer_sum.to_freqs().iter().map(|f| f / depth).collect();

        // M-step: freqs -> copy_nums
        let copy_nums_new = freqs_to_copy_nums(cdbg, &freqs, &copy_nums, false);

        // calc forward probability
        let p_total: Prob = reads
            .par_iter()
            .map(|read| phmm.forward_prob(&param, read))
            .product();

        // log out
        println!(
            "{}\t{:.16}\t{:.16}\t{}\t{:?}\t{:?}",
            i,
            depth,
            p_total.to_log_value(),
            cdbg.to_seqs_string(&copy_nums_new),
            copy_nums_new,
            freqs,
        );

        /*
        // difference check
        if copy_nums_new == copy_nums {
            break;
        } else {
            copy_nums = copy_nums_new;
        }
        */
        copy_nums = copy_nums_new;
    }
}

pub fn true_copy_nums_for_em(
    cdbg: &CompressedDBG,
    reads: &[Vec<u8>],
    param: PHMMParams,
    true_copy_nums: &[u32],
    depth: f64,
) {
    // E-step: copy_nums -> freqs
    let phmm = CDbgPHMM::new(cdbg, true_copy_nums.to_vec());
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
    let freqs: Vec<f64> = layer_sum.to_freqs().iter().map(|f| f / depth).collect();

    // calc forward probability
    let p_total: Prob = reads
        .par_iter()
        .map(|read| phmm.forward_prob(&param, read))
        .product();

    // log out
    println!(
        "{}\t{:.16}\t{:.16}\t{}\t{:?}\t{:?}",
        "true",
        depth,
        p_total.to_log_value(),
        cdbg.to_seqs_string(&true_copy_nums),
        true_copy_nums,
        freqs,
    );
}

/// freqs -> copy_nums function, by fitting with gradient descent and MMWC problem
pub fn freqs_to_copy_nums(
    cdbg: &CompressedDBG,
    freqs: &[f64],
    copy_nums_init: &[u32],
    is_verbose: bool,
) -> Vec<u32> {
    let idg = cdbg.to_indexed_digraph();
    let s = BestFreqState::new(&cdbg, &idg, &freqs, copy_nums_init.to_vec());
    let g = GradientDescent::new(100, is_verbose);
    let mut history = g.run(s);

    // DEBUG
    for (i, h) in history.iter().enumerate() {
        info!("{} {:?}", i, h.copy_nums);
    }

    history.pop().unwrap().copy_nums
}

pub fn freqs_vs_true_copy_nums(cdbg: &CompressedDBG, freqs: &[f64], copy_nums_true: &[u32]) {
    let idg = cdbg.to_indexed_digraph();
    let s = BestFreqState::new(&cdbg, &idg, &freqs, copy_nums_true.to_vec());
    let g = GradientDescent::new(1, true);
    g.run_once(s);
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
