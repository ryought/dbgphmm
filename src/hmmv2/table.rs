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
pub use petgraph::graph::NodeIndex;

/// Struct that stores Forward/Backward algorithm result
/// for the given emissions
///
/// the length of the PHMMResult.tables will be
/// equal to the length of emissions
#[derive(Debug, Clone)]
pub struct PHMMResult<S: Storage<Item = Prob>> {
    pub init_table: PHMMTable<S>,
    pub tables: Vec<PHMMTable<S>>,
}

/// Maximum number of continuous deletion (D -> D transition)
/// allowed in pHMM.
pub const MAX_DEL: usize = 4;

/// Struct for storing Forward/Backward intermediate result
/// for an emission.
///
/// Corresponds to a vector `T[node, type]`
/// `node` is either normal or begin or end node.
/// `type` is either Match, Ins, Del.
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

impl<S: Storage<Item = Prob>> std::fmt::Display for PHMMTable<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "\tMatch\tIns\tDel")?;
        // Begin state
        writeln!(f, "Begin\t{}\t{}", self.mb, self.ib)?;
        // Normal states
        let n_nodes = self.m.len();
        for i in 0..n_nodes {
            let v = NodeIndex::new(i);
            writeln!(f, "{}\t{}\t{}\t{}", i, self.m[v], self.i[v], self.d[v])?;
        }
        // End state
        writeln!(f, "End\t{}", self.e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
}
