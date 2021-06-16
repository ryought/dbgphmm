use super::base::*;
use super::linear::LinearPHMM;
use super::params::PHMMParams;
use crate::prob::Prob;
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;

pub fn test(r: &[u8], q: &[u8]) {
    println!("hmm test for {:?} {:?}", r, q);
    let model = LinearPHMM::from(r);
    for node in model.nodes().iter() {
        println!(
            "{:?} {:?} {:?} {:?} {:?} {}",
            node,
            model.childs(node),
            model.parents(node),
            model.copy_num(node),
            model.emission(node) as char,
            model.init_prob(node),
        );
    }
    let param = PHMMParams::new(
        Prob::from_prob(0.01),
        Prob::from_prob(0.01),
        Prob::from_prob(0.01),
        10,
    );
    let p = model.forward_prob(&param, q);
    println!("prob = {}", p);
}

pub fn test_static() {
    eprintln!("static");
    let model = LinearPHMM::from(b"ATCGATTCGATTAGCT");
    let param = PHMMParams::new(
        Prob::from_prob(0.01),
        Prob::from_prob(0.01),
        Prob::from_prob(0.01),
        10,
    );
    println!("{}", model.as_dot());
    model.sample(&param, 10, 0);
    // let p = model.forward_prob(&param, q);
    // println!("prob = {}", p);
}

pub fn test_random() {
    println!("random testing");
    // let x: u8 = random();
    // println!("{}", x);

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
    let choices = [
        (State::Match, Prob::from_prob(0.9)),
        (State::Ins, Prob::from_prob(0.05)),
        (State::Del, Prob::from_prob(0.05)),
    ];
    let (mut n_a, mut n_b, mut n_c) = (0, 0, 0);
    for _ in 0..100 {
        // let c = &choices.choose_weighted(&mut rng, |item| item.1).unwrap().0;
        // println!("rand {}", rng.gen::<u8>());
        let c = pick_with_prob(&mut rng, &choices);
        println!("choice {:?}", c);
        match c {
            State::Match => n_a += 1,
            State::Ins => n_b += 1,
            State::Del => n_c += 1,
            State::InsBegin => {}
            State::MatchBegin => {}
        }
    }
    println!("choices {} {} {}", n_a, n_b, n_c);
}
