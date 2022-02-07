//!
//! Sampled states and emissions
//!
use super::{Emission, State};
use crate::common::Sequence;

pub fn get_emission_sequence(sampled: &[(State, Emission)]) -> Sequence {
    sampled
        .iter()
        .filter_map(|(_, emission)| match emission {
            Emission::Base(base) => Some(base),
            Emission::Empty => None,
        })
        .copied()
        .collect()
}
