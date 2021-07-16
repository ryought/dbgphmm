use super::annealer::{simple_run, Annealer, SAState};
use super::base::ScoreableState;
use super::grad::GDState;
use crate::compressed_dbg::CompressedDBG;
use crate::cycles::CycleDirection;
use crate::dbg::{DbgHash, DBG};
use crate::hmm;
use crate::hmm::base::PHMM;
use crate::hmm::params::PHMMParams;
use crate::prob::Prob;
use log::{debug, info, warn};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::fmt::Write as FmtWrite;

/// cdbg+phmm with f64 real freqs
#[derive(Clone, PartialEq)]
pub struct FCDbgState<'a> {
    cdbg: &'a CompressedDBG,
    freqs: Vec<f64>,
    reads: Option<&'a [Vec<u8>]>,
    param: PHMMParams,
    score_cache: Option<Prob>,
    parallel: bool,
    delta: f64,
}

impl<'a> FCDbgState<'a> {
    /// create state with given frequencies
    pub fn new(
        cdbg: &'a CompressedDBG,
        freqs: Vec<f64>,
        reads: Option<&'a [Vec<u8>]>,
        param: PHMMParams,
        parallel: bool,
        delta: f64,
    ) -> FCDbgState<'a> {
        FCDbgState {
            cdbg,
            freqs,
            reads,
            param,
            score_cache: None,
            parallel,
            delta,
        }
    }
    /// initial state with all-zero frequencies
    pub fn init(
        cdbg: &'a CompressedDBG,
        reads: Option<&'a [Vec<u8>]>,
        param: PHMMParams,
        parallel: bool,
        delta: f64,
    ) -> FCDbgState<'a> {
        let freqs = vec![0.0; cdbg.n_kmers()];
        FCDbgState::new(cdbg, freqs, reads, param, parallel, delta)
    }
    pub fn forward_score(&self) -> Prob {
        match self.reads {
            Some(reads) => {
                let phmm = hmm::fdbg::FCDbgPHMM::new(self.cdbg, self.freqs.clone());

                if self.parallel {
                    reads
                        .par_iter()
                        .map(|read| phmm.forward_prob(&self.param, read))
                        .product()
                } else {
                    reads
                        .iter()
                        .map(|read| phmm.forward_prob(&self.param, read))
                        .product()
                }
            }
            None => Prob::from_prob(1.0),
        }
    }
    fn calc_score(&mut self) -> Prob {
        match self.score_cache {
            Some(score) => score,
            _ => {
                let score = self.forward_score();
                self.score_cache = Some(score);
                score
            }
        }
    }
}

impl<'a> ScoreableState for FCDbgState<'a> {
    /// Calc the posterior probability
    /// Now it returns only the prior score (no read information)
    fn score(&self) -> f64 {
        self.score_cache.unwrap().to_log_value()
    }
    fn fill_score(&mut self) -> f64 {
        self.calc_score();
        self.score()
    }
    fn as_string(&self) -> String {
        let mut s = String::new();
        match self.score_cache {
            Some(score) => {
                write!(
                    &mut s,
                    "{}\t{}\t{}\t{:?}",
                    0.0,
                    score.to_log_value(),
                    self.cdbg.total_emitable_freq(&self.freqs),
                    self.freqs,
                );
            }
            _ => {
                write!(
                    &mut s,
                    "0\tDefered\t{}\t{:?}",
                    self.cdbg.total_emitable_freq(&self.freqs),
                    self.freqs,
                );
            }
        }
        s
    }
}

impl<'a> GDState for FCDbgState<'a> {
    /// increase/decrease the freq of one k-mer
    fn neighbors(&self) -> Vec<FCDbgState<'a>> {
        let mut neighbors = Vec::new();

        for i in 0..self.cdbg.n_kmers() {
            for is_up in [true, false].iter() {
                let mut new_freqs = self.freqs.to_vec();
                if *is_up {
                    new_freqs[i] += self.delta;
                } else {
                    new_freqs[i] -= self.delta;
                }

                // accept if all freqs are above zero
                if new_freqs.iter().all(|&freq| freq >= 0.0) {
                    let neighbor = FCDbgState {
                        freqs: new_freqs,
                        score_cache: None,
                        param: self.param.clone(),
                        ..*self
                    };
                    neighbors.push(neighbor);
                }
            }
        }

        neighbors
    }
    fn is_duplicate(&self, other: &FCDbgState) -> bool {
        self.freqs == other.freqs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple() {
        let seqs = vec![b"ATTCGATCGATTT".to_vec()];
        let (cdbg, copy_nums) = CompressedDBG::from_seqs(&seqs, 8);

        let freqs: Vec<f64> = copy_nums.iter().map(|&cn| cn as f64).collect();
        let param = PHMMParams::default();
        let mut s = FCDbgState::new(&cdbg, freqs, Some(&seqs), param, true, 0.1);
        let score = s.forward_score();
        println!("score={}", score);
        for n in s.neighbors() {
            println!("{}", n.as_string());
        }

        let param = PHMMParams::default();
        let cycle_vec = cdbg.cycle_vec_from_copy_nums(&copy_nums);
        let mut s_old = crate::optimizer::cdbg::CDbgState::new(
            &cdbg,
            copy_nums,
            cycle_vec,
            10,
            10,
            Some(&seqs),
            param,
            true,
        );
        let score_old = s_old.forward_score();
        println!("score_old={}", score_old);

        assert_eq!(score, score_old);
    }
}
