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
use crate::vector::{NodeVec, Storage};

pub struct PHMMResult<S: Storage<Item = Prob>> {
    pub init_table: PHMMTable<S>,
    pub tables: Vec<PHMMTable<S>>,
}

/// Maximum number of continuous deletion (D -> D transition)
/// allowed in pHMM.
pub const MAX_DEL: usize = 4;

#[derive(Debug, Clone)]
pub struct PHMMTable<S: Storage<Item = Prob>> {
    /// Match node probability
    pub m: NodeVec<S>,
    /// Ins node probability
    pub i: NodeVec<S>,
    /// Del node probability
    pub d: NodeVec<S>,
    /// Match node in begin state
    pub mb: Prob,
    /// Ins node in begin state
    pub ib: Prob,
    /// end state probability
    pub e: Prob,
}

impl<S: Storage<Item = Prob>> PHMMTable<S> {
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

#[cfg(test)]
mod tests {
    use super::*;
}
