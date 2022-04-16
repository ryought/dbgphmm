//!
//! Sampling emissions from the PHMMModel
//!
//! ## Todos
//!
//! * create a special struct for storeing the sampling result instead of `Vec<(State, Emission)>`
//!
use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use crate::common::{Freq, Reads, Sequence};
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;
pub use petgraph::graph::{EdgeIndex, NodeIndex};
use picker::{pick_ins_emission, pick_match_emission, pick_with_prob};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
pub mod history;
pub mod picker;
pub mod utils;
pub use history::{History, Historys};

///
/// HMM hidden states
///
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum State {
    Match(NodeIndex),
    Ins(NodeIndex),
    Del(NodeIndex),
    MatchBegin,
    InsBegin,
    End,
}

impl State {
    ///
    /// Convert State (either Match/Ins/Del(NodeIndex) MatchBegin/InsBegin/End)
    /// into the wrapped NodeIndex.
    /// If state is begin/end, it returns None.
    ///
    pub fn to_node_index(&self) -> Option<NodeIndex> {
        match self {
            State::Match(v) | State::Ins(v) | State::Del(v) => Some(*v),
            _ => None,
        }
    }
}

impl std::fmt::Display for State {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            State::Match(node) => write!(f, "Match({})", node.index()),
            State::Ins(node) => write!(f, "Ins({})", node.index()),
            State::Del(node) => write!(f, "Del({})", node.index()),
            State::MatchBegin => write!(f, "MatchBegin"),
            State::InsBegin => write!(f, "InsBegin"),
            State::End => write!(f, "End"),
        }
    }
}

///
/// HMM emission
///
#[derive(Debug, Copy, Clone)]
pub enum Emission {
    Base(u8),
    Empty,
}

impl Emission {
    ///
    /// Check if the emission is actual base or not.
    ///
    pub fn is_base(&self) -> bool {
        match self {
            Emission::Base(_) => true,
            Emission::Empty => false,
        }
    }
}

impl std::fmt::Display for Emission {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Emission::Base(base) => write!(f, "{}", *base as char),
            Emission::Empty => write!(f, "-"),
        }
    }
}

/// Specs for read sampling
///
#[derive(Clone, Debug)]
pub struct SampleProfile {
    ///
    /// how many reads will be sampled?
    ///
    pub read_amount: ReadAmount,
    ///
    /// seed of rng
    ///
    pub seed: u64,
    ///
    /// length (i.e. bases) of reads
    ///
    pub length: usize,
    ///
    /// start points of reads.
    ///
    pub start_points: StartPoints,
    ///
    /// is endable in the middle of the model.
    ///
    pub endable: bool,
}

///
/// Specifying read amount
///
#[derive(Clone, Debug)]
pub enum ReadAmount {
    /// by the number of reads.
    Count(usize),
    // /// by depth (coverage) of reads.
    // Depth(Freq),
    /// by the total bases of reads.
    TotalBases(usize),
}

///
/// enum for storing start points.
#[derive(Clone, Debug)]
pub enum StartPoints {
    /// all nodes in graph can be a start point.
    /// pick a start node uniformly random
    Random,
    /// all endpoint node in graph
    AllStartPoints,
    /// allow only custom node index
    Custom(Vec<NodeIndex>),
}

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Generate a one-shot sequence of Emission and Hidden states
    /// by running a profile HMM.
    /// Random number generator will be created from the seed.
    ///
    pub fn sample(&self, length: usize, seed: u64) -> History {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        self.sample_rng(&mut rng, length, true)
    }
    ///
    /// Generate a one-shot sequence of Emission and Hidden states
    /// by running a profile HMM.
    /// Random number generator will be created from the seed.
    ///
    pub fn sample_many(&self, length: usize, seed: u64, n_samples: usize) -> Historys {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        Historys(
            (0..n_samples)
                .map(|_| self.sample_rng(&mut rng, length, true))
                .collect(),
        )
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
    ///
    /// Generate reads from sampling from profile
    ///
    pub fn sample_by_profile(&self, profile: &SampleProfile) -> Historys {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(profile.seed);
        let mut sampler = || match &profile.start_points {
            StartPoints::Custom(nodes) => {
                self.sample_rng_from_nodes(&mut rng, profile.length, nodes)
            }
            StartPoints::Random => self.sample_rng(&mut rng, profile.length, true),
            StartPoints::AllStartPoints => {
                panic!("StartPoints::AllStartPoints is not resolved until sample_by_profile")
            }
        };
        let historys = match profile.read_amount {
            ReadAmount::Count(n_reads) => (0..n_reads).map(|_| sampler()).collect(),
            ReadAmount::TotalBases(required_total_bases) => {
                let mut historys = Vec::new();
                let mut total_bases = 0;
                while total_bases < required_total_bases {
                    let history = sampler();
                    total_bases += history.total_bases();
                    historys.push(history);
                }
                historys
            }
        };
        Historys(historys)
    }
    ///
    /// Generate a sequence of Emission and Hidden states
    /// by running a profile HMM using the given rng(random number generator).
    ///
    pub fn sample_rng<R: Rng>(&self, rng: &mut R, length: usize, endable: bool) -> History {
        self.sample_rng_from(rng, length, State::MatchBegin, endable)
    }
    ///
    /// Generate full length sample
    ///
    /// * from
    ///     list of index of the starting nodes
    ///
    pub fn sample_rng_from_nodes<R: Rng>(
        &self,
        rng: &mut R,
        length: usize,
        froms: &[NodeIndex],
    ) -> History {
        let from = self
            .pick_init_node_from(rng, froms)
            .expect("froms are empty");
        self.sample_rng_from(rng, length, State::Match(from), false)
    }
    ///
    /// Generate a sequence of Emission and Hidden states
    /// by running a profile HMM using the given rng(random number generator).
    ///
    /// * from
    ///     a state the sampling starts from. Usually it will set to MatchBegin or the starting node.
    /// * endable
    ///     if true, it allows a transition to End state while sampling, so
    ///     the sampled sequence can be shorter than the given length.
    ///
    pub fn sample_rng_from<R: Rng>(
        &self,
        rng: &mut R,
        length: usize,
        from: State,
        endable: bool,
    ) -> History {
        let mut history = History::new();

        // (1) init
        let mut state = from;
        let mut emission;
        match from {
            State::MatchBegin => {}
            _ => {
                emission = self.make_emission(rng, from);
                history.push(state, emission);
            }
        }

        for _ in 0..length {
            // (2) step
            (state, emission) = self.sample_step(rng, state, endable);
            history.push(state, emission);

            if let State::End = state {
                break;
            }
        }

        history
    }
    fn sample_step<R: Rng>(&self, rng: &mut R, now: State, endable: bool) -> (State, Emission) {
        // (1) transition to the next state
        let next = self.make_transition(rng, now, endable);
        // (2) emission from the next state
        let emission = self.make_emission(rng, next);

        (next, emission)
    }
    ///
    /// do a transition from `now: State` to `next: State` in PHMM
    ///
    /// * now
    /// * endable
    ///
    fn make_transition<R: Rng>(&self, rng: &mut R, now: State, endable: bool) -> State {
        let param = &self.param;
        let p_end = if endable {
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
                            (State::End, p_end),
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
                            (State::End, p_end),
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
    ///
    /// randomly-picks a transition to a initial node
    /// from the choices given.
    ///
    fn pick_init_node_from<R: Rng>(&self, rng: &mut R, nodes: &[NodeIndex]) -> Option<NodeIndex> {
        let choices: Vec<(NodeIndex, Prob)> = nodes
            .iter()
            .map(|&v| (v, self.node(v).init_prob()))
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
    use crate::common::{ni, sequence_to_string};
    use crate::hmmv2::mocks::mock_linear_phmm;

    #[test]
    fn hmm_sample_state() {
        let s = State::Match(NodeIndex::new(10));
        println!("{:?}", s);
    }
    #[test]
    fn hmm_sample_state_to_node_index() {
        assert_eq!(State::Match(ni(10)).to_node_index(), Some(ni(10)));
        assert_eq!(State::Ins(ni(2)).to_node_index(), Some(ni(2)));
        assert_eq!(State::InsBegin.to_node_index(), None);
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
        println!("{}", hist);
        assert_eq!(read, b"CGATC");
        println!("{:?}", sequence_to_string(&read));
    }
    #[test]
    fn hmm_sample_mock_linear_high_error() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(3);
        let phmm = mock_linear_phmm(PHMMParams::high_error());
        let hist = phmm.sample(10, 0);
        let read = hist.to_sequence();
        println!("{}", hist);
        assert_eq!(read, b"CACAACGT");
    }
    #[test]
    fn hmm_sample_mock_linear_full_length() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(3);
        let phmm = mock_linear_phmm(PHMMParams::high_error());
        println!("{}", phmm);
        for i in 0..10 {
            println!("{}", i);
            // start node is ni(0)
            let hist = phmm.sample_rng_from_nodes(&mut rng, 100, &[NodeIndex::new(0)]);
            println!("{}", hist);
            // first is ni(0)
            assert_eq!(hist.0[0].0.to_node_index(), Some(NodeIndex::new(0)));
            // 2nd-last is ni(9)
            assert_eq!(
                hist.0[hist.0.len() - 2].0.to_node_index(),
                Some(NodeIndex::new(9))
            );
        }
    }
    #[test]
    fn hmm_sample_mock_linear_picker_from() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
        let phmm = mock_linear_phmm(PHMMParams::default());
        println!("{}", phmm);

        // if choices are only node2
        for i in 0..10 {
            let node = phmm.pick_init_node_from(&mut rng, &[ni(2)]);
            assert_eq!(node, Some(ni(2)));
        }

        // if choices are only node0 and node1
        let mut count = (0, 0);
        for i in 0..100 {
            let node = phmm.pick_init_node_from(&mut rng, &[ni(0), ni(1)]);
            assert!(node.is_some());
            let v = node.unwrap();
            assert!(v == ni(0) || v == ni(1));
            println!("{:?}", v);
            if v == ni(0) {
                count.0 += 1;
            } else if v == ni(1) {
                count.1 += 1;
            };
        }
        println!("{:?}", count);
        assert!(40 <= count.0 && count.0 <= 60);
        assert!(40 <= count.1 && count.1 <= 60);
    }
    #[test]
    fn hmm_sample_mock_linear_by_profile() {
        let phmm = mock_linear_phmm(PHMMParams::default());

        // a
        let start_point = ni(3);
        let hists = phmm.sample_by_profile(&SampleProfile {
            read_amount: ReadAmount::Count(10),
            seed: 0,
            length: 100,
            start_points: StartPoints::Custom(vec![start_point]),
            endable: false,
        });
        for hist in hists.iter() {
            println!("a {} {}", hist, hist.to_sequence().len());
            assert_eq!(hist.0[0].0, State::Match(start_point));
            assert_eq!(hist.to_sequence().len(), 7);
        }

        // b
        let hists = phmm.sample_by_profile(&SampleProfile {
            read_amount: ReadAmount::TotalBases(500),
            seed: 0,
            length: 100,
            start_points: StartPoints::Random,
            endable: false,
        });
        for hist in hists.iter() {
            println!("b {}", sequence_to_string(&hist.to_sequence()));
        }
        let total_bases: usize = hists.iter().map(|hist| hist.total_bases()).sum();
        println!("total={}", total_bases);
        assert_eq!(total_bases, 500);
    }
}
