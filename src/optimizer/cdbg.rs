use super::base::{simple_run, SAState};
use crate::compressed_dbg::CompressedDBG;
use crate::dbg::{DbgHash, DBG};
use crate::distribution::normal_bin;
use log::{debug, info, warn};
use rand::prelude::*;
use std::fmt::Write as FmtWrite;

/// SAState for cdbg
/// for cycle debugging purpose
#[derive(Clone)]
struct CDbgState<'a> {
    cdbg: &'a CompressedDBG,
    copy_nums: Vec<u32>,
    ave_size: u32,
    std_size: u32,
}

impl<'a> CDbgState<'a> {
    /// create state with given copy-nums
    fn new(cdbg: &CompressedDBG, copy_nums: Vec<u32>, ave_size: u32, std_size: u32) -> CDbgState {
        CDbgState {
            cdbg,
            copy_nums,
            ave_size,
            std_size,
        }
    }
    /// initial state with all-zero copy-nums
    fn init(cdbg: &CompressedDBG, ave_size: u32, std_size: u32) -> CDbgState {
        let copy_nums = vec![0; cdbg.n_kmers()];
        CDbgState {
            cdbg,
            copy_nums,
            ave_size,
            std_size,
        }
    }
}

impl<'a> SAState for CDbgState<'a> {
    /// Calc the posterior probability
    /// Now it returns only the prior score (no read information)
    fn score(&self) -> f64 {
        normal_bin(
            self.cdbg.total_emitable_copy_num(&self.copy_nums),
            self.ave_size,
            self.std_size,
        )
        .to_log_value()
    }
    /// Pick cycles randomly and return new state
    fn next<R: Rng>(&self, rng: &mut R) -> CDbgState<'a> {
        // 1. pick a cycle
        let cycle_id = rng.gen_range(0..self.cdbg.n_cycles());
        debug!("next: {}", cycle_id);
        // 2. pick a direction (+1 or -1) of modifying the copy-nums
        let min_copy_num = self
            .cdbg
            .cycle_components(cycle_id)
            .iter()
            .map(|v| self.copy_nums[v.0])
            .min()
            .unwrap();
        let is_down: bool = min_copy_num > 0 && rng.gen();
        debug!("direction: {} (min: {})", is_down, min_copy_num);
        let mut copy_nums = self.copy_nums.clone();
        for v in self.cdbg.cycle_components(cycle_id) {
            if is_down {
                copy_nums[v.0] -= 1;
            } else {
                copy_nums[v.0] += 1;
            }
        }
        CDbgState {
            cdbg: self.cdbg,
            copy_nums,
            ave_size: self.ave_size,
            std_size: self.std_size,
        }
    }
    fn as_string(&self) -> String {
        let mut s = String::new();
        writeln!(
            &mut s,
            "{:?} {} {}",
            self.copy_nums,
            self.cdbg.is_consistent_copy_num(&self.copy_nums),
            self.cdbg.total_emitable_copy_num(&self.copy_nums)
        );
        s
    }
}

pub fn test() {
    let mut d = DbgHash::new();
    d.add_seq(b"ATCGATTCGATCGATTCGATAGATCG", 8);
    let cdbg = CompressedDBG::from(&d, 8);

    let s0 = CDbgState::init(&cdbg, 26, 10);
    simple_run(s0, 100);
}
