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
use super::grad::GDState;
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
    #[allow(unused_must_use)]
    fn as_string(&self) -> String {
        let mut s = String::new();
        write!(
            &mut s,
            "{}\t{}\t{}\t{}\t{:?}",
            self.score(),
            0.0,
            self.cdbg.total_emitable_copy_num(&self.copy_nums),
            self.cdbg.to_seqs_string(&self.copy_nums),
            self.copy_nums
        );
        s
    }
}

impl<'a> SAState for FreqState<'a> {
    /// Suggest next move by these two strategies
    /// - unbiased move
    /// - biased move
    fn next<R: Rng>(&self, rng: &mut R) -> FreqState<'a> {
        /*
        for (cycle_id, is_up) in self
            .cdbg
            .cycle_and_direction_candidates(&self.copy_nums)
            .into_iter()
        {
            println!("{}\t{}", self.score(), s.as_string());
        }
        */
        let candidates_with_score: Vec<((usize, bool), f64)> = self
            .cdbg
            .cycle_and_direction_candidates(&self.copy_nums)
            .into_iter()
            .map(|(cycle_id, is_up)| {
                let copy_nums = self.cdbg.update_by_cycle(&self.copy_nums, cycle_id, is_up);
                let cycle_vec =
                    self.cdbg
                        .update_cycle_vec_by_cycle(&self.cycle_vec, cycle_id, is_up);
                let s = FreqState {
                    cdbg: self.cdbg,
                    freqs: self.freqs,
                    copy_nums,
                    cycle_vec,
                };
                ((cycle_id, is_up), s.score())
            })
            .collect();
        for x in candidates_with_score.iter() {
            println!("{} {} {}", x.0 .0, x.0 .1, x.1);
        }

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

impl<'a> GDState for FreqState<'a> {
    fn neighbors(&self) -> Vec<FreqState<'a>> {
        let mut neighbors = Vec::new();

        for cycle_id in 0..self.cdbg.n_cycles() {
            for is_up in [true, false].iter() {
                if self.cdbg.is_acceptable(&self.copy_nums, cycle_id, *is_up) {
                    let copy_nums = self.cdbg.update_by_cycle(&self.copy_nums, cycle_id, *is_up);
                    let cycle_vec =
                        self.cdbg
                            .update_cycle_vec_by_cycle(&self.cycle_vec, cycle_id, *is_up);

                    let neighbor = FreqState {
                        cdbg: self.cdbg,
                        freqs: self.freqs,
                        copy_nums,
                        cycle_vec,
                    };
                    neighbors.push(neighbor);
                }
            }
        }

        neighbors
    }
    fn is_duplicate(&self, other: &FreqState) -> bool {
        self.copy_nums == other.copy_nums
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
