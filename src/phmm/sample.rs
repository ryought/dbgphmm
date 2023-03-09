//!
//! Sampling emissions from PHMM
//!

// struct to store the sampled states and emissions
pub mod history;
pub use history::History;

// pick an element from list
pub mod picker;
use picker::{pick_ins_emission, pick_match_emission, pick_with_prob};

use super::{Emission, State, PHMM};

use crate::prob::Prob;

pub use petgraph::graph::{EdgeIndex, NodeIndex};
use rand::prelude::*;

///
/// Starting condition
///
/// * Random
/// * Custom
///
#[derive(Clone, Debug)]
pub enum StartPoints {
    /// all nodes in graph can be a start point.
    /// pick a start node uniformly random
    Random,
    /// allow only custom node index
    Custom(Vec<NodeIndex>),
}

///
/// End condition
///
/// * NoChild
/// * StateCount
/// * VisitEnd
///
#[derive(Clone, Debug, Copy)]
pub enum EndCond {
    ///
    /// Never stop transition, except there is no child.
    ///
    NoChild,
    ///
    /// By number of state transition
    ///
    StateCount(usize),
    ///
    /// By number of emission
    ///
    EmitCount(usize),
    ///
    /// By transition to End node by probability of p_end.
    ///
    VisitEnd,
}

impl EndCond {
    pub fn is_endable(&self) -> bool {
        match self {
            EndCond::VisitEnd => true,
            _ => false,
        }
    }
}

///
/// Public functions
///
impl PHMM {
    ///
    ///
    ///
    pub fn sample_rng<R: Rng>(&self, rng: &mut R, start: &StartPoints, end: &EndCond) -> History {
        self.sample_rng_from(rng, start, end, State::MatchBegin)
    }
    ///
    /// Generate a sequence of Emission and Hidden states
    /// by running a profile HMM using the given rng(random number generator).
    ///
    pub fn sample_rng_from<R: Rng>(
        &self,
        rng: &mut R,
        start: &StartPoints,
        end: &EndCond,
        from: State,
    ) -> History {
        let mut history = History::new();

        // counter
        let mut n_state = 0;
        let mut n_emission = 0;

        // (1) init
        let mut state = from;
        let mut emission;
        match from {
            State::MatchBegin => {}
            _ => {
                emission = self.make_emission(rng, from);
                history.push(state, emission);

                // update counter
                n_emission += 1;
            }
        }

        loop {
            // (2) step
            (state, emission) = self.sample_step(rng, state, start, end);
            history.push(state, emission);

            // update counter
            n_state += 1;
            if emission.is_base() {
                n_emission += 1;
            }

            // loop exit conditions
            if let State::End = state {
                break;
            }
            match end {
                EndCond::StateCount(c) => {
                    if n_state == *c {
                        break;
                    }
                }
                EndCond::EmitCount(c) => {
                    if n_emission == *c {
                        break;
                    }
                }
                EndCond::NoChild => continue,
                EndCond::VisitEnd => continue,
            }
        }

        history
    }
}

//
// internal functions
//
// * make Transition and Emission
//
impl PHMM {
    ///
    /// Make Transition and Emission from current state `now: State` in PHMM
    ///
    fn sample_step<R: Rng>(
        &self,
        rng: &mut R,
        now: State,
        start: &StartPoints,
        end: &EndCond,
    ) -> (State, Emission) {
        // (1) transition to the next state
        let next = self.make_transition(rng, now, start, end);
        // (2) emission from the next state
        let emission = self.make_emission(rng, next);

        (next, emission)
    }
    ///
    /// Transition from `now: State` to `next: State` in PHMM
    ///
    fn make_transition<R: Rng>(
        &self,
        rng: &mut R,
        now: State,
        start: &StartPoints,
        end: &EndCond,
    ) -> State {
        let param = &self.param;
        let p_end = if end.is_endable() {
            param.p_end
        } else {
            Prob::from_prob(0.0)
        };
        match now {
            State::Match(node) => {
                let child = self.pick_child(rng, node);
                match child {
                    Some(child) => {
                        let choices = [
                            (State::Match(child), param.p_MM),
                            (State::Ins(node), param.p_MI),
                            (State::Del(child), param.p_MD),
                            (State::End, p_end),
                        ];
                        pick_with_prob(rng, &choices).unwrap()
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
                            (State::End, p_end),
                        ];
                        pick_with_prob(rng, &choices).unwrap()
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
                            (State::End, p_end),
                        ];
                        pick_with_prob(rng, &choices).unwrap()
                    }
                    None => State::End,
                }
            }
            State::MatchBegin => {
                let node = self.pick_init_node(rng, start);
                match node {
                    Some(node) => {
                        let choices = [
                            (State::InsBegin, param.p_MI),
                            (State::Match(node), param.p_MM),
                            (State::Del(node), param.p_MD),
                        ];
                        pick_with_prob(rng, &choices).unwrap()
                    }
                    None => State::End,
                }
            }
            State::InsBegin => {
                let node = self.pick_init_node(rng, start);
                match node {
                    Some(node) => {
                        let choices = [
                            (State::InsBegin, param.p_II),
                            (State::Match(node), param.p_IM),
                            (State::Del(node), param.p_ID),
                        ];
                        pick_with_prob(rng, &choices).unwrap()
                    }
                    None => State::End,
                }
            }
            State::End => unreachable!(),
        }
    }
    ///
    /// Emission from the current state `state: State` in PHMM
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
    /// Pick an initial node in the graph
    ///
    fn pick_init_node<R: Rng>(&self, rng: &mut R, start: &StartPoints) -> Option<NodeIndex> {
        let choices: Vec<(NodeIndex, Prob)> = match start {
            StartPoints::Random => {
                // all nodes
                self.nodes()
                    .map(|(v, vw)| (v, vw.init_prob()))
                    .filter(|(_, prob)| !prob.is_zero())
                    .collect()
            }
            StartPoints::Custom(nodes) => {
                // the specified nodes only
                nodes
                    .iter()
                    .map(|&v| (v, self.init_prob(v)))
                    .filter(|(_, prob)| !prob.is_zero())
                    .collect()
            }
        };

        pick_with_prob(rng, &choices)
    }
    ///
    /// Pick a child node from current node
    ///
    fn pick_child<R: Rng>(&self, rng: &mut R, node: NodeIndex) -> Option<NodeIndex> {
        let choices: Vec<(NodeIndex, Prob)> = self
            .childs(node)
            .map(|(_, child, ew)| (child, ew.trans_prob()))
            .filter(|(_, prob)| !prob.is_zero())
            .collect();

        pick_with_prob(rng, &choices)
    }
}
