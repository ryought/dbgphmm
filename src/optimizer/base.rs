//! simulated annealing optimization module
//!

use log::{debug, info, warn};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::fmt::Write as FmtWrite;

pub trait SAState: Clone {
    /// Score of the state. Typically, this score corresponds to
    /// the log likelihood (log p) of the model to be maximized.
    fn score(&self) -> f64;
    /// get a randomly-picked neighbor state (using rng)
    fn next<R: Rng>(&self, rng: &mut R) -> Self;
    /// output for logging
    fn as_string(&self) -> String;
}

/// find SAState with global maximum score
pub struct Annealer {
    init_temp: f64,
    cooling_rate: f64,
}

impl Annealer {
    pub fn new() -> Annealer {
        Annealer {
            init_temp: 1.0,
            cooling_rate: 0.8,
        }
    }
    pub fn temp_schedule(&self, iteration: u64) -> f64 {
        self.init_temp * self.cooling_rate.powi(iteration as i32)
    }
    pub fn step<T: SAState, R: Rng>(
        &self,
        rng: &mut R,
        now: &T,
        iteration: u64,
        is_verbose: bool,
    ) -> T {
        let temp = self.temp_schedule(iteration);
        let next = now.next(rng);
        // TODO avoid recalculation of now.score()
        let p = ((next.score() - now.score()) / temp).exp();

        // accept next state with prob min(1, p)
        let is_accepted = rng.gen_bool(p.min(1f64));
        if is_verbose {
            println!(
                "{}\t{:.32}\t{:.32}\t{:.32}\t{:.32}\t{}\t{}",
                iteration,
                temp,
                now.score(),
                next.score(),
                p,
                is_accepted,
                now.as_string(),
            )
        }
        if is_accepted {
            next
        } else {
            now.clone()
        }
    }
    pub fn run<T: SAState, R: Rng>(&self, rng: &mut R, init_state: T, n_iteration: u64) -> Vec<T> {
        let mut states = Vec::new();
        let mut now = init_state;
        for i in 0..n_iteration {
            let new_state = self.step(rng, &now, i, false);
            states.push(now);
            now = new_state;
        }
        states
    }
    pub fn run_with_log<T: SAState, R: Rng>(
        &self,
        rng: &mut R,
        init_state: T,
        n_iteration: u64,
    ) -> Vec<T> {
        let mut states = Vec::new();
        let mut now = init_state;
        for i in 0..n_iteration {
            let new_state = self.step(rng, &now, i, true);
            states.push(now);
            now = new_state;
        }
        states
    }
}

/// greedy transition of SAStates for debugging
pub fn simple_run<T: SAState>(init: T, n_iteration: u64) {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
    let mut now = init;
    for i in 0..n_iteration {
        debug!("now: score={} state={}", now.score(), now.as_string());
        now = now.next(&mut rng);
    }
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

impl SAState for TestState {
    /// score by how near the target value
    fn score(&self) -> f64 {
        -(0.3 - (self.bits.iter().sum::<u8>() as f64 / self.bits.len() as f64)).abs()
    }
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
    fn as_string(&self) -> String {
        let mut s = String::new();
        writeln!(&mut s, "{:?}", self.bits);
        s
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
    let a = Annealer::new();
    let history = a.run(&mut rng, init, 100);
    for state in history.iter() {
        println!("f: {:?}", state.bits);
    }
}
