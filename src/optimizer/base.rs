//! score maximization module

use super::annealer::{Annealer, SAState};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::fmt::Write as FmtWrite;

pub trait ScoreableState {
    /// Score of the state. Typically, this score corresponds to
    /// the log likelihood (log p) of the model to be maximized.
    fn score(&self) -> f64;
    /// fill the score on cache
    fn fill_score(&mut self) -> f64;
    /// output for logging
    fn as_string(&self) -> String;
}

/// random vector
/// with norm regulation
#[derive(Debug, Clone)]
struct TestState {
    bits: Vec<u8>,
}

impl TestState {
    fn new(n: usize) -> TestState {
        TestState { bits: vec![0; n] }
    }
}

impl ScoreableState for TestState {
    /// score by how near the target value
    fn score(&self) -> f64 {
        -(0.3 - (self.bits.iter().sum::<u8>() as f64 / self.bits.len() as f64)).abs()
    }
    fn fill_score(&mut self) -> f64 {
        self.score()
    }
    fn as_string(&self) -> String {
        let mut s = String::new();
        writeln!(&mut s, "{:?}", self.bits);
        s
    }
}

impl SAState for TestState {
    /// pick random position and swap it
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

pub fn test() {
    let mut s = TestState::new(30);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
    for i in 0..100 {
        let score = s.score();
        println!("#{} state={:?}, score={}", i, s.bits, score);
        s = s.next(&mut rng);
    }

    let init = TestState::new(30);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
    let a = Annealer::new(1.0, 0.8);
    let history = a.run(&mut rng, init, 100);
    for state in history.iter() {
        println!("f: {:?}", state.bits);
    }
}
