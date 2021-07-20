//! Frequency only optimizer
//!
//! It is (cdbg, freqs) -> (copy_nums)
//! i.e. find the best copy_nums that is closest to freqs
//! by minimizing differences between freqs and copy_nums.
//!
//! possible solutions for this problem are
//!
//! - **on integer space** uses simulated annealing
//! - **on real space** relaxation and stregnthen the effect of constraints

use super::annealer::{Annealer, SAState};
use super::base::ScoreableState;
use crate::compressed_dbg::CompressedDBG;
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::fmt::Write as FmtWrite;

#[derive(Clone, PartialEq)]
pub struct FreqState<'a> {
    cdbg: &'a CompressedDBG,
    freqs: &'a [f64],
    copy_nums: Vec<u32>,
    cycle_vec: Vec<u32>,
}

impl<'a> FreqState<'a> {
    pub fn new(cdbg: &'a CompressedDBG, freqs: &'a [f64], copy_nums: &'a [u32]) -> FreqState<'a> {
        let cycle_vec = cdbg.cycle_vec_from_copy_nums(copy_nums);
        FreqState {
            cdbg,
            freqs,
            copy_nums: copy_nums.to_vec(),
            cycle_vec,
        }
    }
    pub fn init(cdbg: &'a CompressedDBG, freqs: &'a [f64]) -> FreqState<'a> {
        let copy_nums = vec![0; cdbg.n_kmers()];
        let cycle_vec = vec![0; cdbg.n_cycles()];
        FreqState {
            cdbg,
            freqs,
            copy_nums,
            cycle_vec,
        }
    }
    fn choose_cycle_and_direction<R: Rng>(&self, rng: &mut R) -> (usize, bool) {
        self.cdbg
            .cycle_and_direction_candidates(&self.copy_nums)
            .choose(rng)
            .unwrap()
            .clone()
    }
}

impl<'a> ScoreableState for FreqState<'a> {
    /// calc score of the state
    fn score(&self) -> f64 {
        let diff: f64 = self
            .copy_nums
            .iter()
            .zip(self.freqs.iter())
            .map(|(&cn, &f)| (cn as f64 - f).powi(2))
            .sum();
        diff * -1f64
    }
    fn fill_score(&mut self) -> f64 {
        // not using a cache, so simply returns score()
        self.score()
    }
    fn as_string(&self) -> String {
        let mut s = String::new();
        write!(&mut s, "{}\t{:?}", self.score(), self.copy_nums);
        s
    }
}

impl<'a> SAState for FreqState<'a> {
    /// Suggest next move by these two strategies
    /// - unbiased move
    /// - biased move
    fn next<R: Rng>(&self, rng: &mut R) -> FreqState<'a> {
        let (cycle_id, is_up) = self.choose_cycle_and_direction(rng);
        let copy_nums = self.cdbg.update_by_cycle(&self.copy_nums, cycle_id, is_up);
        let cycle_vec = self
            .cdbg
            .update_cycle_vec_by_cycle(&self.cycle_vec, cycle_id, is_up);
        FreqState {
            cdbg: self.cdbg,
            freqs: self.freqs,
            copy_nums,
            cycle_vec,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_freq_state() {
        let seqs = vec![b"ATTCGATCGATTT".to_vec()];
        let (cdbg, copy_nums) = CompressedDBG::from_seqs(&seqs, 8);
        let freqs = cdbg.copy_nums_to_freqs(&copy_nums);
        let s = FreqState::init(&cdbg, &freqs);
        println!("{}", s.as_string());

        let a = Annealer::new(1.0, 0.8);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
        a.run_with_log(&mut rng, s, 100);
    }
}
