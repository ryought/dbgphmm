//!
//! Markov chain monte carlo methods
//!
//! Random sample
//!

trait Model {
    fn score(&self) -> f64;
}

struct Runner {
    // proposal: P,
}
