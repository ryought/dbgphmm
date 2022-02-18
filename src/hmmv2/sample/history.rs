//!
//! Struct for storeing sampling results
//!
use super::super::freq::NodeFreqs;
use super::super::trans_table::EdgeFreqs;
use super::{Emission, State};
use crate::common::Sequence;

///
/// Struct for storing sampling results from HMM
///
pub struct History(Vec<(State, Emission)>);

impl History {
    ///
    /// Constructor of empty sample store.
    ///
    pub fn new() -> Self {
        History(Vec::new())
    }
    ///
    /// Append a new state and its emission
    ///
    pub fn push(&mut self, state: State, emission: Emission) {
        self.0.push((state, emission));
    }
    ///
    /// Create base sequence `Vec<u8>` from sampling history
    ///
    pub fn to_sequence(&self) -> Sequence {
        self.0
            .iter()
            .filter_map(|(_, emission)| match emission {
                Emission::Base(base) => Some(base),
                Emission::Empty => None,
            })
            .copied()
            .collect()
    }
    /// TODO
    ///
    pub fn to_node_freqs(&self) -> NodeFreqs {
        unimplemented!();
    }
    /// TODO
    ///
    pub fn to_edge_freqs(&self) -> EdgeFreqs {
        unimplemented!();
    }
}

impl std::fmt::Display for History {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for (state, emission) in self.0.iter() {
            writeln!(f, "{} -> {}", state, emission)?;
        }
        Ok(())
    }
}
