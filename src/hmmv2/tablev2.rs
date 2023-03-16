//!
//! PHMMTable V2
//!
//! Updated version of PHMMTable using SparseVec
//!
use crate::prob::{p, Prob};
use petgraph::graph::NodeIndex;
use sparsevec::SparseVec;

///
/// Maximum number of continuous deletion (D -> D transition)
/// allowed in pHMM.
///
pub const MAX_DEL: usize = 4;

///
/// Maximum number of active nodes (= nodes which have high probability)
///
pub const MAX_ACTIVE_NODES: usize = 400;

/// NodeVec is `SparseVec<Prob, NodeIndex, SIZE>`
///
/// `NodeVec x[node] = prob`
///
pub type NodeVec<const N: usize> = SparseVec<Prob, NodeIndex, N>;

///
/// PHMMTable V2
///
/// A Vector `T[Node, Type]` to compute forward/backward DP.
/// * Node is normal graph node and End node
/// * Type is Match/Ins/Del
///
/// Implementation:
/// * m[node], i[node], d[node]: SparseVec<Node, Prob>
/// * mb, ib, e: Prob
///
#[derive(Debug, Clone)]
pub struct PHMMTable<const N: usize> {
    ///
    /// Match node probability
    /// `Table[Match, node]`
    ///
    pub m: NodeVec<N>,
    ///
    /// Ins node probability
    /// `Table[Ins, node]`
    ///
    pub i: NodeVec<N>,
    ///
    /// Del node probability
    /// `Table[Del, node]`
    ///
    pub d: NodeVec<N>,
    ///
    /// Match node in begin state
    /// `Table[Match, Begin]`
    ///
    pub mb: Prob,
    ///
    /// Ins node in begin state
    /// `Table[Ins, Begin]`
    ///
    pub ib: Prob,
    ///
    /// end state probability
    /// `Table[End]`
    ///
    pub e: Prob,
}

///
/// PHMMTables
///
/// Vec of PHMMTable to store tables of emissions
///
#[derive(Debug, Clone)]
pub struct PHMMTables<const N: usize> {
    pub init_table: PHMMTable<N>,
    pub tables: Vec<PHMMTable<N>>,
    pub is_forward: bool,
}
