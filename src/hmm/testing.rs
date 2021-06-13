use super::base::{PHMMParams, PHMM};
use super::linear::LinearPHMM;
use crate::prob::Prob;

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
