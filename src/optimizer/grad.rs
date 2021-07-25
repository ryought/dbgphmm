//! gradient-based greedy optimization module
//!
//! Start from the seed state

use super::base::ScoreableState;
use log::{debug, info, warn};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::fmt::Write as FmtWrite;

pub trait GDState: Clone + ScoreableState {
    fn neighbors(&self) -> Vec<Self>;
    fn is_duplicate(&self, other: &Self) -> bool;
}

/// find local maximum score by gradient descent
/// 1. start from seed
/// 2. calc all scores of neighbors of current state
/// 3-a. pick the most highest state
/// 3-b. if all scores are lower than current state, stop iteration
pub struct GradientDescent {
    max_iteration: u64,
    is_verbose: bool,
    stop_when_loop: bool,
}

impl GradientDescent {
    pub fn new(max_iteration: u64, is_verbose: bool) -> GradientDescent {
        GradientDescent {
            max_iteration,
            is_verbose,
            stop_when_loop: true,
        }
    }
    /// one step of climbing
    /// calc scores of all neighbors
    pub fn step<T: GDState>(&self, now: &T, iteration: u64) -> Option<T> {
        let neighbor_scores: Vec<(T, f64)> = now
            .neighbors()
            .into_iter()
            .enumerate()
            .map(|(i, mut neighbor)| {
                neighbor.fill_score();
                let score = neighbor.score();
                // log output
                if self.is_verbose {
                    println!(
                        "{}\t{}\t{:.32}\t{:.32}\t{:.32}\t{}\t{}\t{}",
                        iteration,
                        i,
                        now.score(),
                        neighbor.score(),
                        0.0,
                        false,
                        now.as_string(),
                        neighbor.as_string(),
                    )
                }
                (neighbor, score)
            })
            .collect();
        if neighbor_scores.len() == 0 {
            None
        } else {
            let max_neighbor = neighbor_scores
                .into_iter()
                .filter(|(_, x)| !x.is_nan())
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap());
            match max_neighbor {
                Some((neighbor, score)) => {
                    // new score is better than current score
                    if score > now.score() {
                        Some(neighbor)
                    } else {
                        None
                    }
                }
                None => None,
            }
        }
    }
    pub fn run<T: GDState>(&self, init_state: T) -> Vec<T> {
        let mut states = Vec::new();
        let mut now = init_state;
        now.fill_score();
        states.push(now.clone());
        if self.is_verbose {
            println!("# INIT\t{:.32}\t{}", now.score(), now.as_string());
        }
        for iteration in 0..self.max_iteration {
            match self.step(&now, iteration) {
                Some(next) => {
                    // loop check
                    if self.stop_when_loop && states.iter().any(|s| next.is_duplicate(s)) {
                        if self.is_verbose {
                            println!("# LOOP\t{:.32}\t{}", next.score(), next.as_string());
                        }
                        break;
                    }
                    if self.is_verbose {
                        println!("# MOVE\t{:.32}\t{}", next.score(), next.as_string());
                    }
                    states.push(now);
                    now = next;
                }
                None => {
                    if self.is_verbose {
                        println!("# STOP\t{:.32}\t{}", now.score(), now.as_string());
                    }
                    break;
                }
            }
        }
        states
    }
    pub fn run_once<T: GDState>(&self, mut state: T) -> T {
        state.fill_score();
        println!("# ONCE\t{:.32}\t{}", state.score(), state.as_string());
        state
    }
}
