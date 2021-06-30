use crate::cycles::DbgTree;
use crate::dbg::{DbgHash, DBG};

/// SAState for dbg
#[derive(Debug, Clone)]
struct DbgCycleState {
    dbg: &DbgHash,
    tree: &DbgTree,
}

impl DbgCycleState {
    fn new(dbg: &DbgPHMM) -> DbgCycleState {
        DbgCycleState { bits: vec![0; n] }
    }
}

impl SAState for DbgCycleState {
    /// calc the prior probability of the genome size
    fn score(&self) -> f64 {}
    /// pick cycles randomly and modify copy-nums accordingly
    fn next<R: Rng>(&self, rng: &mut R) -> TestState {
        let mut new_bits = self.bits.clone();
        let x = new_bits.choose_mut(rng).unwrap();
        if *x > 0 {
            *x = 0;
        } else {
            *x = 1;
        }
        TestState { bits: new_bits }
    }
}
