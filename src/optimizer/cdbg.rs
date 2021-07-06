use super::base::{simple_run, Annealer, SAState};
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
use std::fmt::Write as FmtWrite;

/// SAState for cdbg and cdbg-based phmm
#[derive(Clone)]
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
        }
    }
    /// initial state with all-zero copy-nums
    pub fn init(
        cdbg: &'a CompressedDBG,
        ave_size: u32,
        std_size: u32,
        reads: Option<&'a [Vec<u8>]>,
        param: PHMMParams,
    ) -> CDbgState<'a> {
        let copy_nums = vec![0; cdbg.n_kmers()];
        let cycle_vec = vec![0; cdbg.n_cycles()];
        CDbgState::new(cdbg, copy_nums, cycle_vec, ave_size, std_size, reads, param)
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
    fn forward_score(&self) -> Prob {
        // TODO omit this cloning of copy_nums(Cdbgphmm does not have to real vec)
        match self.reads {
            Some(reads) => {
                let phmm = hmm::cdbg::CDbgPHMM::new(self.cdbg, self.copy_nums.clone());

                reads
                    .iter()
                    .enumerate()
                    .map(|(i, read)| {
                        // warn!("read i={}", i);
                        phmm.forward_prob(&self.param, read)
                    })
                    .product()
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

impl<'a> SAState for CDbgState<'a> {
    /// Calc the posterior probability
    /// Now it returns only the prior score (no read information)
    fn score(&self) -> f64 {
        (self.prior_score_cache.unwrap() * self.forward_score_cache.unwrap()).to_log_value()
    }
    fn fill_score(&mut self) -> f64 {
        self.calc_score();
        self.score()
    }
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
        }
    }
    fn as_string(&self) -> String {
        let mut s = String::new();
        match (self.prior_score_cache, self.forward_score_cache) {
            (Some(prior_score), Some(forward_score)) => {
                write!(
                    &mut s,
                    "{}\t{}\t{}\t{:?}",
                    prior_score,
                    forward_score,
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

pub fn test() {
    let mut d = DbgHash::new();
    d.add_seq(b"ATCGATTCGATCGATTCGATAGATCG", 8);
    let cdbg = CompressedDBG::from(&d, 8);

    let param = PHMMParams::default();
    let init = CDbgState::init(&cdbg, 26, 10, None, param);
    // simple_run(init, 100);

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
    let a = Annealer::new(1.0, 0.8);
    let history = a.run(&mut rng, init, 100);
}
