use super::base::{Node, PHMM};
use super::params::PHMMParams;
use crate::prob::Prob;
use log::info;
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;

#[derive(Debug, Copy, Clone)]
pub enum State {
    Match,
    Ins,
    Del,
    MatchBegin,
    InsBegin,
}

pub trait PHMMSampler: PHMM {
    // sampling related
    fn sample(&self, param: &PHMMParams, length: u32, seed: u64) -> Vec<u8> {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        let mut emissions: Vec<u8> = Vec::new();
        let mut now: (State, Node) = (State::MatchBegin, Node(0));
        for i in 0..=length {
            info!("iter {} {:?}", i, now);
            // 1. emission
            match now {
                (State::Match, v) => {
                    let emission_match_choices = emission_match_choices(&param, self.emission(&v));
                    let emission = pick_with_prob(&mut rng, &emission_match_choices);
                    info!("emit {} <- {}", emission as char, self.emission(&v) as char);
                    emissions.push(emission);
                }
                (State::Ins, _) | (State::InsBegin, _) => {
                    let emission_ins_choices = emission_ins_choices(&param);
                    let emission = pick_with_prob(&mut rng, &emission_ins_choices);
                    info!("emit {} <- 0", emission as char);
                    emissions.push(emission);
                }
                _ => {}
            };
            // 2. transition state and kmer
            let next = match now {
                (State::Match, v) => {
                    let state_choices = [
                        (State::Match, param.p_MM),
                        (State::Ins, param.p_MI),
                        (State::Del, param.p_MD),
                    ];
                    let state = pick_with_prob(&mut rng, &state_choices);
                    match state {
                        State::Ins => (state, v),
                        _ => {
                            let move_choices: Vec<(Node, Prob)> = self
                                .childs(&v)
                                .into_iter()
                                .map(|w| (w, self.trans_prob(&v, &w)))
                                .collect();
                            let w = pick_with_prob(&mut rng, &move_choices);
                            (state, w)
                        }
                    }
                }
                (State::Ins, v) => {
                    let state_choices = [
                        (State::Match, param.p_IM),
                        (State::Ins, param.p_II),
                        (State::Del, param.p_ID),
                    ];
                    let state = pick_with_prob(&mut rng, &state_choices);
                    match state {
                        State::Ins => (state, v),
                        _ => {
                            let move_choices: Vec<(Node, Prob)> = self
                                .childs(&v)
                                .into_iter()
                                .map(|w| (w, self.trans_prob(&v, &w)))
                                .collect();
                            let w = pick_with_prob(&mut rng, &move_choices);
                            (state, w)
                        }
                    }
                }
                (State::Del, v) => {
                    let state_choices = [
                        (State::Match, param.p_DM),
                        (State::Ins, param.p_DI),
                        (State::Del, param.p_DD),
                    ];
                    let state = pick_with_prob(&mut rng, &state_choices);
                    match state {
                        State::Ins => (state, v),
                        _ => {
                            let move_choices: Vec<(Node, Prob)> = self
                                .childs(&v)
                                .into_iter()
                                .map(|w| (w, self.trans_prob(&v, &w)))
                                .collect();
                            let w = pick_with_prob(&mut rng, &move_choices);
                            (state, w)
                        }
                    }
                }
                (State::MatchBegin, _) => {
                    let state_choices = [
                        (State::InsBegin, param.p_MI),
                        (State::Match, param.p_MM),
                        (State::Del, param.p_MD),
                    ];
                    let state = pick_with_prob(&mut rng, &state_choices);
                    match state {
                        State::InsBegin => (state, Node(0)),
                        _ => {
                            let move_choices: Vec<(Node, Prob)> = self
                                .nodes()
                                .into_iter()
                                .map(|w| (w, self.init_prob(&w)))
                                .collect();
                            let w = pick_with_prob(&mut rng, &move_choices);
                            (state, w)
                        }
                    }
                }
                (State::InsBegin, _) => {
                    let state_choices = [
                        (State::InsBegin, param.p_II),
                        (State::Match, param.p_IM),
                        (State::Del, param.p_ID),
                    ];
                    let state = pick_with_prob(&mut rng, &state_choices);
                    match state {
                        State::InsBegin => (state, Node(0)),
                        _ => {
                            let move_choices: Vec<(Node, Prob)> = self
                                .nodes()
                                .into_iter()
                                .map(|w| (w, self.init_prob(&w)))
                                .collect();
                            let w = pick_with_prob(&mut rng, &move_choices);
                            (state, w)
                        }
                    }
                }
            };
            now = next;
        }
        emissions
    }
}

pub fn pick_with_prob<R: Rng, T: Copy>(rng: &mut R, choices: &[(T, Prob)]) -> T {
    choices
        .choose_weighted(rng, |item| item.1.to_value())
        .unwrap()
        .0
}

pub fn emission_ins_choices(param: &PHMMParams) -> Vec<(u8, Prob)> {
    let choices = [b'A', b'C', b'G', b'T']
        .iter()
        .map(|&base| (base, param.p_random))
        .collect();
    choices
}

pub fn emission_match_choices(param: &PHMMParams, emission: u8) -> Vec<(u8, Prob)> {
    let choices = [b'A', b'C', b'G', b'T']
        .iter()
        .map(|&base| {
            (
                base,
                if base == emission {
                    param.p_match
                } else {
                    param.p_mismatch
                },
            )
        })
        .collect();
    choices
}
