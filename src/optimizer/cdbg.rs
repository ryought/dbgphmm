use super::base::{simple_run, Annealer, SAState};
use crate::compressed_dbg::CompressedDBG;
use crate::cycles::CycleDirection;
use crate::dbg::{DbgHash, DBG};
use log::{debug, info, warn};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::fmt::Write as FmtWrite;

/// SAState for cdbg
/// for cycle debugging purpose
#[derive(Clone)]
pub struct CDbgState<'a> {
    cdbg: &'a CompressedDBG,
    copy_nums: Vec<u32>,
    cycle_vec: Vec<u32>, // how many times cycle basis was used?
    ave_size: u32,
    std_size: u32,
}

impl<'a> CDbgState<'a> {
    /// create state with given copy-nums
    pub fn new(
        cdbg: &CompressedDBG,
        copy_nums: Vec<u32>,
        cycle_vec: Vec<u32>,
        ave_size: u32,
        std_size: u32,
    ) -> CDbgState {
        CDbgState {
            cdbg,
            copy_nums,
            cycle_vec,
            ave_size,
            std_size,
        }
    }
    /// initial state with all-zero copy-nums
    pub fn init(cdbg: &CompressedDBG, ave_size: u32, std_size: u32) -> CDbgState {
        let copy_nums = vec![0; cdbg.n_kmers()];
        let cycle_vec = vec![0; cdbg.n_cycles()];
        CDbgState {
            cdbg,
            copy_nums,
            cycle_vec,
            ave_size,
            std_size,
        }
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
}

impl<'a> SAState for CDbgState<'a> {
    /// Calc the posterior probability
    /// Now it returns only the prior score (no read information)
    fn score(&self) -> f64 {
        self.cdbg
            .prior_score(&self.copy_nums, self.ave_size, self.std_size)
            .to_log_value()
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
        }
    }
    fn as_string(&self) -> String {
        let mut s = String::new();
        write!(
            &mut s,
            "{}\t{:?}",
            self.cdbg.total_emitable_copy_num(&self.copy_nums),
            self.cycle_vec,
        );
        s
    }
}

pub fn test() {
    let mut d = DbgHash::new();
    d.add_seq(b"ATCGATTCGATCGATTCGATAGATCG", 8);
    let cdbg = CompressedDBG::from(&d, 8);

    let init = CDbgState::init(&cdbg, 26, 10);
    // simple_run(init, 100);

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
    let a = Annealer::new();
    let history = a.run(&mut rng, init, 100);
    for state in history.iter() {
        println!("f: {:?}", state.score());
    }
}
