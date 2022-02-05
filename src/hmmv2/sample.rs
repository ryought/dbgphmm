//!
//! Sampling emissions from the PHMMModel
//!
use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use crate::hmm::params::PHMMParams;
use crate::prob::Prob;
pub use petgraph::graph::{EdgeIndex, NodeIndex};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
pub mod picker;

#[derive(Debug, Copy, Clone)]
pub enum State {
    Match(NodeIndex),
    Ins(NodeIndex),
    Del(NodeIndex),
    MatchBegin,
    InsBegin,
    End,
}

#[derive(Debug, Copy, Clone)]
pub struct Emission(Option<u8>);

/// valid bases
const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    pub fn sample(&self, length: usize) {}
    fn sample_init<R: Rng>(&self, rng: &mut R) -> (State, Emission) {
        (State::MatchBegin, Emission(None))
    }
    fn sample_step<R: Rng>(&self, rng: &mut R, now: State) -> (State, Emission) {
        let param = &self.param;
        // (1) transition
        let state = match now {
            State::Match(node) => {
                let child = self.pick_child(rng, node);
                match child {
                    Some(child) => {
                        let choices = [
                            (State::Match(child), param.p_MM),
                            (State::Ins(node), param.p_MI),
                            (State::Del(child), param.p_MD),
                            (State::End, param.p_end),
                        ];
                        pick_with_prob(rng, &choices)
                    }
                    None => State::End,
                }
            }
            State::Ins(node) => {
                let child = self.pick_child(rng, node);
                match child {
                    Some(child) => {
                        let choices = [
                            (State::Match(child), param.p_IM),
                            (State::Ins(node), param.p_II),
                            (State::Del(child), param.p_ID),
                            (State::End, param.p_end),
                        ];
                        pick_with_prob(rng, &choices)
                    }
                    None => State::End,
                }
            }
            State::Del(node) => {
                let child = self.pick_child(rng, node);
                match child {
                    Some(child) => {
                        let choices = [
                            (State::Match(child), param.p_DM),
                            (State::Ins(node), param.p_DI),
                            (State::Del(child), param.p_DD),
                            (State::End, param.p_end),
                        ];
                        pick_with_prob(rng, &choices)
                    }
                    None => State::End,
                }
            }
            State::MatchBegin => {
                let node = self.pick_init_node(rng);
                match node {
                    Some(node) => {
                        let choices = [
                            (State::InsBegin, param.p_MI),
                            (State::Match(node), param.p_MM),
                            (State::Del(node), param.p_MD),
                        ];
                        pick_with_prob(rng, &choices)
                    }
                    None => State::End,
                }
            }
            State::InsBegin => {
                let node = self.pick_init_node(rng);
                match node {
                    Some(node) => {
                        let choices = [
                            (State::InsBegin, param.p_II),
                            (State::Match(node), param.p_IM),
                            (State::Del(node), param.p_ID),
                        ];
                        pick_with_prob(rng, &choices)
                    }
                    None => State::End,
                }
            }
            State::End => panic!(),
        };

        // (2) emission
        let emission = match state {
            State::Match(node) => pick_match_emission(rng, self.emission(node), param),
            State::Ins(_) | State::InsBegin => pick_ins_emission(rng, param),
            _ => Emission(None),
        };

        (state, emission)
    }
    fn pick_child<R: Rng>(&self, rng: &mut R, node: NodeIndex) -> Option<NodeIndex> {
        let choices: Vec<(NodeIndex, Prob)> = self
            .childs(node)
            .map(|(_, child, ew)| (child, ew.trans_prob()))
            .filter(|(_, prob)| !prob.is_zero())
            .collect();

        if choices.len() == 0 {
            // there is no child (with >0 transition probability)
            None
        } else {
            Some(pick_with_prob(rng, &choices))
        }
    }
    fn pick_init_node<R: Rng>(&self, rng: &mut R) -> Option<NodeIndex> {
        // TODO
        // when start from head mode, this should be go to the first node
        let choices: Vec<(NodeIndex, Prob)> = self
            .nodes()
            .map(|(v, vw)| (v, vw.init_prob()))
            .filter(|(_, prob)| !prob.is_zero())
            .collect();

        if choices.len() == 0 {
            // there is no child (with >0 transition probability)
            None
        } else {
            Some(pick_with_prob(rng, &choices))
        }
    }
}

///
/// pick randomly from the choices with its own probability.
///
pub fn pick_with_prob<R: Rng, T: Copy>(rng: &mut R, choices: &[(T, Prob)]) -> T {
    choices
        .choose_weighted(rng, |item| item.1.to_value())
        .unwrap()
        .0
}

///
/// Pick an emission at `Ins` state, according to the PHMM params.
/// It selects uniformaly random on `ACGT`
///
pub fn pick_ins_emission<R: Rng>(rng: &mut R, param: &PHMMParams) -> Emission {
    let base = *BASES
        .choose_weighted(rng, |_base| param.p_random.to_value())
        .unwrap();
    Emission(Some(base))
}

///
/// Pick an emission at `Match` state (with a state emission), according to the PHMM params.
/// It selects a state emission with probability `p_match`, and other emissions with probability
/// `p_mismatch`.
///
pub fn pick_match_emission<R: Rng>(rng: &mut R, emission: u8, param: &PHMMParams) -> Emission {
    let base = *BASES
        .choose_weighted(rng, |&base| {
            if base == emission {
                param.p_match.to_value()
            } else {
                param.p_mismatch.to_value()
            }
        })
        .unwrap();
    Emission(Some(base))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hmm_sample_state() {
        let s = State::Match(NodeIndex::new(10));
        println!("{:?}", s);
    }
}
