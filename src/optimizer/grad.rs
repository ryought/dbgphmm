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
}

/// find local maximum score by gradient descent
/// 1. start from seed
/// 2. calc all scores of neighbors of current state
/// 3-a. pick the most highest state
/// 3-b. if all scores are lower than current state, stop iteration
pub struct GradientDescent {
    max_iteration: u64,
    is_verbose: bool,
}

impl GradientDescent {
    pub fn new(max_iteration: u64, is_verbose: bool) -> GradientDescent {
        GradientDescent {
            max_iteration,
            is_verbose,
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
                Some((neighbor, _)) => Some(neighbor),
                None => None,
            }
        }
    }
    pub fn run<T: GDState>(&self, init_state: T) -> Vec<T> {
        let mut states = Vec::new();
        let mut now = init_state;
        now.fill_score();
        for iteration in 0..self.max_iteration {
            match self.step(&now, iteration) {
                Some(next) => {
                    if self.is_verbose {
                        println!("# MOVE\t{:.32}\t{}", next.score(), next.as_string());
                    }
                    states.push(now);
                    now = next;
                }
                None => {
                    println!("# STOP");
                }
            }
        }
        states
    }
}
