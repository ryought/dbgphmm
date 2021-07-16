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

/// SAState for cdbg and cdbg-based phmm
#[derive(Clone, PartialEq)]
pub struct CDbgState<'a> {
    cdbg: &'a CompressedDBG,
    pub copy_nums: Vec<u32>,
    cycle_vec: Vec<u32>, // how many times cycle basis was used?
    ave_size: u32,
    std_size: u32,
    reads: Option<&'a [Vec<u8>]>,
    param: PHMMParams,
    prior_score_cache: Option<Prob>,
    forward_score_cache: Option<Prob>,
    parallel: bool,
}

impl<'a> CDbgState<'a> {
    /// create state with given copy-nums
    pub fn new(
        cdbg: &'a CompressedDBG,
        copy_nums: Vec<u32>,
        cycle_vec: Vec<u32>,
        ave_size: u32,
        std_size: u32,
        reads: Option<&'a [Vec<u8>]>,
        param: PHMMParams,
        parallel: bool,
    ) -> CDbgState<'a> {
        CDbgState {
            cdbg,
            copy_nums,
            cycle_vec,
            ave_size,
            std_size,
            reads,
            param,
            prior_score_cache: None,
            forward_score_cache: None,
            parallel,
        }
    }
    /// initial state with all-zero copy-nums
    pub fn init(
        cdbg: &'a CompressedDBG,
        ave_size: u32,
        std_size: u32,
        reads: Option<&'a [Vec<u8>]>,
        param: PHMMParams,
        parall: bool,
    ) -> CDbgState<'a> {
        let copy_nums = vec![0; cdbg.n_kmers()];
        let cycle_vec = vec![0; cdbg.n_cycles()];
        CDbgState::new(
            cdbg, copy_nums, cycle_vec, ave_size, std_size, reads, param, parall,
        )
    }
    /// generate states randomly by adding random cycle basis
    pub fn random<R: Rng>(
        rng: &mut R,
        n_basis: u32,
        cdbg: &'a CompressedDBG,
        ave_size: u32,
        std_size: u32,
        reads: Option<&'a [Vec<u8>]>,
        param: PHMMParams,
        parall: bool,
    ) -> CDbgState<'a> {
        let mut copy_nums = vec![0; cdbg.n_kmers()];
        let mut cycle_vec = vec![0; cdbg.n_cycles()];
        for _ in 0..n_basis {
            let (cycle_id, is_up) = cdbg
                .cycle_and_direction_candidates(&copy_nums)
                .choose(rng)
                .unwrap()
                .clone();
            copy_nums = cdbg.update_by_cycle(&copy_nums, cycle_id, is_up);
            cycle_vec = cdbg.update_cycle_vec_by_cycle(&cycle_vec, cycle_id, is_up);
        }
        CDbgState::new(
            cdbg, copy_nums, cycle_vec, ave_size, std_size, reads, param, parall,
        )
    }
    fn choose_cycle_and_direction<R: Rng>(&self, rng: &mut R) -> (usize, bool) {
        for _ in 0..100 {
            // 1. pick a cycle
            let cycle_id = rng.gen_range(0..self.cdbg.n_cycles());

            // 2. pick a direction (+1 or -1) of modifying the copy-nums
            let is_up_movable = self.cdbg.is_acceptable(&self.copy_nums, cycle_id, true);
            let is_down_movable = self.cdbg.is_acceptable(&self.copy_nums, cycle_id, false);
            match (is_up_movable, is_down_movable) {
                (true, true) => return (cycle_id, rng.gen()),
                (true, false) => return (cycle_id, true),
                (false, true) => return (cycle_id, false),
                _ => continue,
            };
        }
        panic!("movable cycle not found");
    }
    fn prior_score(&self) -> Prob {
        self.cdbg
            .prior_score(&self.copy_nums, self.ave_size, self.std_size)
    }
    pub fn forward_score(&self) -> Prob {
        // TODO omit this cloning of copy_nums(Cdbgphmm does not have to real vec)
        match self.reads {
            Some(reads) => {
                let phmm = hmm::cdbg::CDbgPHMM::new(self.cdbg, self.copy_nums.clone());

                if self.parallel {
                    reads
                        .par_iter()
                        .map(|read| phmm.forward_prob(&self.param, read))
                        .product()
                } else {
                    reads
                        .iter()
                        .enumerate()
                        .map(|(i, read)| {
                            // warn!("read i={}", i);
                            phmm.forward_prob(&self.param, read)
                        })
                        .product()
                }
            }
            None => Prob::from_prob(1.0),
        }
    }
    fn calc_score(&mut self) -> Prob {
        match (self.prior_score_cache, self.forward_score_cache) {
            (Some(prior_score), Some(forward_score)) => prior_score * forward_score,
            _ => {
                let prior_score = self.prior_score();
                let forward_score = self.forward_score();
                self.prior_score_cache = Some(prior_score);
                self.forward_score_cache = Some(forward_score);
                warn!("score filled {} {}", prior_score, forward_score);
                prior_score * forward_score
            }
        }
    }
}

impl<'a> ScoreableState for CDbgState<'a> {
    /// Calc the posterior probability
    /// Now it returns only the prior score (no read information)
    fn score(&self) -> f64 {
        (self.prior_score_cache.unwrap() * self.forward_score_cache.unwrap()).to_log_value()
    }
    fn fill_score(&mut self) -> f64 {
        self.calc_score();
        self.score()
    }
    fn as_string(&self) -> String {
        let mut s = String::new();
        match (self.prior_score_cache, self.forward_score_cache) {
            (Some(prior_score), Some(forward_score)) => {
                write!(
                    &mut s,
                    "{}\t{}\t{}\t{:?}",
                    prior_score.to_log_value(),
                    forward_score.to_log_value(),
                    self.cdbg.total_emitable_copy_num(&self.copy_nums),
                    self.cycle_vec,
                );
            }
            _ => {
                write!(
                    &mut s,
                    "Defered\tDefered\t{}\t{:?}",
                    self.cdbg.total_emitable_copy_num(&self.copy_nums),
                    self.cycle_vec,
                );
            }
        }
        s
    }
}

impl<'a> SAState for CDbgState<'a> {
    /// Pick cycles randomly and return new state
    fn next<R: Rng>(&self, rng: &mut R) -> CDbgState<'a> {
        let (cycle_id, is_up) = self.choose_cycle_and_direction(rng);
        let copy_nums = self.cdbg.update_by_cycle(&self.copy_nums, cycle_id, is_up);
        let cycle_vec = self
            .cdbg
            .update_cycle_vec_by_cycle(&self.cycle_vec, cycle_id, is_up);
        let is_consistent = self.cdbg.is_consistent_copy_num(&self.copy_nums);
        assert!(is_consistent);
        info!(
            "next cycle={} up={} is_consistent={}",
            cycle_id, is_up, is_consistent
        );
        CDbgState {
            cdbg: self.cdbg,
            copy_nums,
            cycle_vec,
            ave_size: self.ave_size,
            std_size: self.std_size,
            reads: self.reads,
            param: self.param.clone(),
            prior_score_cache: None,
            forward_score_cache: None,
            parallel: self.parallel,
        }
    }
}

impl<'a> GDState for CDbgState<'a> {
    /// change for each basis into two (up/down) directions
    fn neighbors(&self) -> Vec<CDbgState<'a>> {
        let mut neighbors = Vec::new();

        for cycle_id in 0..self.cdbg.n_cycles() {
            for is_up in [true, false].iter() {
                if self.cdbg.is_acceptable(&self.copy_nums, cycle_id, *is_up) {
                    let copy_nums = self.cdbg.update_by_cycle(&self.copy_nums, cycle_id, *is_up);
                    let cycle_vec =
                        self.cdbg
                            .update_cycle_vec_by_cycle(&self.cycle_vec, cycle_id, *is_up);

                    let neighbor = CDbgState {
                        copy_nums,
                        cycle_vec,
                        prior_score_cache: None,
                        forward_score_cache: None,
                        param: self.param.clone(),
                        ..*self
                    };
                    neighbors.push(neighbor);
                }
            }
        }

        neighbors
    }
    fn is_duplicate(&self, other: &CDbgState) -> bool {
        self.copy_nums == other.copy_nums
    }
}

pub fn test() {
    let mut d = DbgHash::new();
    d.add_seq(b"ATCGATTCGATCGATTCGATAGATCG", 8);
    let cdbg = CompressedDBG::from(&d, 8);

    let param = PHMMParams::default();
    let init = CDbgState::init(&cdbg, 26, 10, None, param, false);
    // simple_run(init, 100);

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
    let a = Annealer::new(1.0, 0.8);
    let history = a.run(&mut rng, init, 100);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn random_state() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(10);
        let seqs = vec![b"ATTCGATCGATTT".to_vec()];
        let (cdbg, copy_nums) = CompressedDBG::from_seqs(&seqs, 8);
        let param = PHMMParams::default();

        let s0 = CDbgState::random(&mut rng, 10, &cdbg, 10, 10, None, param, false);
        assert!(cdbg.is_consistent_copy_num(&s0.copy_nums));
    }
}
