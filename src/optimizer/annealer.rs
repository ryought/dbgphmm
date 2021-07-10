//! simulated annealing optimization module
//!

use super::base::ScoreableState;
use log::{debug, info, warn};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::fmt::Write as FmtWrite;

pub trait SAState: Clone + ScoreableState {
    /// get a randomly-picked neighbor state (using rng)
    fn next<R: Rng>(&self, rng: &mut R) -> Self;
}

/// find SAState with global maximum score
pub struct Annealer {
    init_temp: f64,
    cooling_rate: f64,
}

impl Annealer {
    pub fn new(init_temp: f64, cooling_rate: f64) -> Annealer {
        Annealer {
            init_temp,
            cooling_rate,
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
        let mut next = now.next(rng);
        next.fill_score();
        // TODO avoid recalculation of now.score()
        let p = ((next.score() - now.score()) / temp).exp();

        // accept next state with prob min(1, p)
        let is_accepted = rng.gen_bool(p.min(1f64));
        if is_verbose {
            println!(
                "{}\t{:.32}\t{:.32}\t{:.32}\t{:.32}\t{}\t{}\t{}",
                iteration,
                temp,
                now.score(),
                next.score(),
                p,
                is_accepted,
                now.as_string(),
                next.as_string(),
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
        now.fill_score();
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
        now.fill_score();
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
