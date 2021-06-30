/// simulated annealing optimization module
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;

trait SAState: Clone + std::fmt::Debug {
    /// score of the state
    fn score(&self) -> f64;
    /// get a randomly-picked neighbor state (using rng)
    fn next<R: Rng>(&self, rng: &mut R) -> Self;
}

/// find SAState with global minimum score
struct Annealer {
    init_temp: f64,
    cooling_rate: f64,
    // rng: &mut R,
    // now: T,
    // history: Vec<T>,
}

impl Annealer {
    fn new() -> Annealer {
        Annealer {
            init_temp: 1.0,
            cooling_rate: 0.8,
        }
    }
    fn temp_schedule(&self, iteration: u64) -> f64 {
        self.init_temp * self.cooling_rate.powi(iteration as i32)
    }
    fn step<T: SAState, R: Rng>(&self, rng: &mut R, now: &T, iteration: u64) -> T {
        let temp = self.temp_schedule(iteration);
        let next = now.next(rng);
        // TODO avoid recalculation of now.score()
        let p = ((next.score() - now.score()) / temp).exp();
        // accept next state with prob min(1, p)
        println!(
            "{} prob={} P(next)={} P(now)={}",
            iteration,
            p,
            next.score(),
            now.score()
        );
        if rng.gen_bool(p.min(1f64)) {
            println!("A {:?}", next);
            next
        } else {
            println!("R {:?}", now);
            now.clone()
        }
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
