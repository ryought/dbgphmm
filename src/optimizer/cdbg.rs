use super::base::{simple_run, SAState};
use crate::compressed_dbg::CompressedDBG;
use crate::dbg::{DbgHash, DBG};
use log::{debug, info, warn};
use rand::prelude::*;
use std::fmt::Write as FmtWrite;

/// SAState for cdbg
/// for cycle debugging purpose
#[derive(Clone)]
struct CDbgState<'a> {
    cdbg: &'a CompressedDBG,
    copy_nums: Vec<u32>,
}

impl<'a> CDbgState<'a> {
    /// create state with given copy-nums
    fn new(cdbg: &CompressedDBG, copy_nums: Vec<u32>) -> CDbgState {
        CDbgState { cdbg, copy_nums }
    }
    /// initial state with all-zero copy-nums
    fn init(cdbg: &CompressedDBG) -> CDbgState {
        let copy_nums = vec![0; cdbg.n_kmers()];
        CDbgState { cdbg, copy_nums }
    }
}

impl<'a> SAState for CDbgState<'a> {
    /// Calc the posterior probability
    /// Now it returns only the prior score (no read information)
    fn score(&self) -> f64 {
        // TODO temporary constant scores
        1.0
    }
    /// Pick cycles randomly and return new state
    fn next<R: Rng>(&self, rng: &mut R) -> CDbgState<'a> {
        let cycle_id = rng.gen_range(0..self.cdbg.n_cycles());
        debug!("next: {}", cycle_id);
        let mut copy_nums = self.copy_nums.clone();
        // TODO how much increase/decrease copy-nums?
        for v in self.cdbg.cycle_components(cycle_id) {
            copy_nums[v.0] += 1;
        }
        CDbgState {
            cdbg: self.cdbg,
            copy_nums,
        }
    }
    fn as_string(&self) -> String {
        let mut s = String::new();
        writeln!(
            &mut s,
            "{:?} {}",
            self.copy_nums,
            self.cdbg.is_consistent_copy_num(&self.copy_nums)
        );
        s
    }
}

pub fn test() {
    let mut d = DbgHash::new();
    d.add_seq(b"ATCGATTCGATCGATTCGATAGATCG", 8);
    let cdbg = CompressedDBG::from(&d, 8);

    let s0 = CDbgState::init(&cdbg);
    simple_run(s0, 100);
}
