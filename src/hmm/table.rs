//!
//! Table definitions
//!
//! ## PHMMTable
//!
//! the prob assigned for each (nodes, node_type)
//!
//! F[Match,v] or B[Match,v]
//!
use super::veclikewrap::NodeVec;
use crate::prob::Prob;
use crate::veclike::{DenseVec, VecLike};

// TODO separate the initial layer?
pub struct PHMMResult<V: VecLike<Prob>>(pub Vec<PHMMTable<V>>);

pub const MAX_DEL: usize = 4;

#[derive(Debug, Clone)]
pub struct PHMMTable<V: VecLike<Prob>> {
    /// Match node probability
    pub m: NodeVec<V>,
    /// Ins node probability
    pub i: NodeVec<V>,
    /// Del node probability
    pub d: NodeVec<V>,
    /// Match node in begin state
    pub mb: Prob,
    /// Ins node in begin state
    pub ib: Prob,
    /// end state probability
    pub e: Prob,
}

impl<V: VecLike<Prob>> PHMMTable<V> {
    pub fn new(n_nodes: usize, m: Prob, i: Prob, d: Prob, mb: Prob, ib: Prob, e: Prob) -> Self {
        PHMMTable {
            m: NodeVec::new(n_nodes, m),
            i: NodeVec::new(n_nodes, i),
            d: NodeVec::new(n_nodes, d),
            mb,
            ib,
            e,
        }
    }
}

// impl PHMMForward for PHMM {}

#[cfg(test)]
mod tests {
    use super::*;
}
