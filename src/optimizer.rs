/// simulated annealing optimization module
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;

trait SAState {
    /// score of the state
    fn score(&self) -> f64;
    /// get a randomly-picked neighbor state (using rng)
    fn next<R: Rng>(&self, rng: &mut R) -> Self;
}

struct Annealer {
    init_temp: f64,
    // rng: &mut R,
    // now: T,
    // history: Vec<T>,
}

impl Annealer {
    fn new() -> Annealer {
        Annealer { init_temp: 10.0 }
    }
    fn temp_schedule(&self, iteration: u64) -> f64 {
        self.init_temp
    }
    fn step<T: SAState, R: Rng>(&self, rng: &mut R, now: &T, iteration: u64) -> T {
        let temp = self.temp_schedule(iteration);
        // accept or reject next state?
        let next = now.next(rng);
        next
    }
    fn run<T: SAState, R: Rng>(&self, rng: &mut R, init_state: T, n_iteration: u64) -> Vec<T> {
        let mut states = Vec::new();
        let mut now = init_state;
        for i in 0..n_iteration {
            let new_state = self.step(rng, &now, i);
            states.push(now);
            now = new_state;
        }
        states
    }
}

/// random vector
/// with norm regulation
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
        0.3 - (self.bits.iter().sum::<u8>() as f64 / self.bits.len() as f64)
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
    let history = a.run(&mut rng, init, 10);
    for state in history.iter() {
        println!("f: {:?}", state.bits);
    }
}
