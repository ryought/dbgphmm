//!
//! Sampling emissions from the PHMMModel
//!
use super::common::{PHMMEdge, PHMMModel, PHMMNode};
pub mod picker;

#[derive(Debug, Copy, Clone)]
pub enum State {
    Match,
    Ins,
    Del,
    MatchBegin,
    InsBegin,
    End,
}

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    pub fn sample(&self, length: usize) {}
}
