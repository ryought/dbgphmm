//!
//! Sampling emissions from the PHMMModel
//!
//! ## Todos
//!
//! * create a special struct for storeing the sampling result instead of `Vec<(State, Emission)>`
//!
use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use crate::common::Sequence;
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;
pub use petgraph::graph::{EdgeIndex, NodeIndex};
use picker::{pick_ins_emission, pick_match_emission, pick_with_prob};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
pub mod history;
pub mod picker;
pub mod utils;
pub use history::History;

///
/// HMM hidden states
///
#[derive(Debug, Copy, Clone)]
pub enum State {
    Match(NodeIndex),
    Ins(NodeIndex),
    Del(NodeIndex),
    MatchBegin,
    InsBegin,
    End,
}

///
/// HMM emission
///
#[derive(Debug, Copy, Clone)]
pub enum Emission {
    Base(u8),
    Empty,
}

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Generate a sequence of Emission and Hidden states
    /// by running a profile HMM using rng(random number generator).
    ///
    pub fn sample(&self, length: usize, seed: u64) -> History {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        let mut history = History::new();

        // (1) init
        let (mut state, _) = self.sample_init(&mut rng);
        let mut emission;

        for _ in 0..length {
            // (2) step
            (state, emission) = self.sample_step(&mut rng, state);
            history.push(state, emission);

            if let State::End = state {
                break;
            }
        }

        history
    }
    ///
    /// Generate a sequence of emissions
    /// by running a profile HMM using rng(random number generator).
    ///
    /// (sampling reads from the model)
    ///
    pub fn sample_read(&self, length: usize, seed: u64) -> Sequence {
        let history = self.sample(length, seed);
        history.to_sequence()
    }
    fn sample_init<R: Rng>(&self, _rng: &mut R) -> (State, Emission) {
        (State::MatchBegin, Emission::Empty)
    }
    fn sample_step<R: Rng>(&self, rng: &mut R, now: State) -> (State, Emission) {
        // (1) transition to the next state
        let next = self.make_transition(rng, now);
        // (2) emission from the next state
        let emission = self.make_emission(rng, next);

        (next, emission)
    }
    ///
    /// do a transition from `now: State` to `next: State` in PHMM
    ///
    fn make_transition<R: Rng>(&self, rng: &mut R, now: State) -> State {
        let param = &self.param;
        match now {
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
        }
    }
    ///
    /// do an emission from the current state `state: State` in PHMM.
    ///
    fn make_emission<R: Rng>(&self, rng: &mut R, state: State) -> Emission {
        let param = &self.param;
        match state {
            State::Match(node) => pick_match_emission(rng, self.emission(node), param),
            State::Ins(_) | State::InsBegin => pick_ins_emission(rng, param),
            _ => Emission::Empty,
        }
    }
    ///
    /// randomly-picks a transition to a child from the current node, according to
    /// the assigned transition probability
    ///
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
    ///
    /// randomly-picks a transition to a initial node, according to
    /// the assigned initial probability from Begin state.
    ///
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::sequence_to_string;
    use crate::hmmv2::mocks::mock_linear_phmm;

    #[test]
    fn hmm_sample_state() {
        let s = State::Match(NodeIndex::new(10));
        println!("{:?}", s);
    }

    #[test]
    fn hmm_sample_mock_linear_picker() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
        let phmm = mock_linear_phmm(PHMMParams::default());
        println!("{}", phmm);

        // pick_init_node
        for _ in 0..100 {
            let node = phmm.pick_init_node(&mut rng);
        }

        for _ in 0..10 {
            // head node
            assert_eq!(
                phmm.pick_child(&mut rng, NodeIndex::new(0)),
                Some(NodeIndex::new(1)),
            );
            // tail node
            assert_eq!(phmm.pick_child(&mut rng, NodeIndex::new(9)), None);
        }
    }
    #[test]
    fn hmm_sample_mock_linear() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(3);
        let phmm = mock_linear_phmm(PHMMParams::zero_error());
        let hist = phmm.sample(5, 0);
        let read = hist.to_sequence();
        println!("{:?}", hist);
        assert_eq!(read, b"CGATC");
        println!("{:?}", sequence_to_string(&read));
    }
}
