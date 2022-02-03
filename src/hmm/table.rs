//!
//! Table definitions
//!
//! ## PHMMTable
//!
//! the prob assigned for each (nodes, node_type)
//!
//! F[Match,v] or B[Match,v]
//!
use crate::prob::Prob;
use crate::veclike::{DenseVec, VecLike};

// TODO separate the initial layer?
pub struct PHMMResult<V: VecLike<Prob>>(pub Vec<PHMMTable<V>>);

pub const MAX_DEL: usize = 4;

#[derive(Debug, Clone)]
pub struct PHMMTable<V: VecLike<Prob>> {
    /// Match node probability
    pub m: V,
    /// Ins node probability
    pub i: V,
    /// Del node probability
    pub d: V,
    /// Match node in begin state
    pub mb: Prob,
    /// Ins node in begin state
    pub ib: Prob,
    /// end state probability
    pub e: Prob,
}

// impl PHMMForward for PHMM {}
