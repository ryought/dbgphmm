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
pub type NodeVec<const M: usize> = SparseVec<Prob, NodeIndex, M>;

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
pub struct PHMMTable<const M: usize> {
    ///
    /// Match node probability
    /// `Table[Match, node]`
    ///
    pub m: NodeVec<M>,
    ///
    /// Ins node probability
    /// `Table[Ins, node]`
    ///
    pub i: NodeVec<M>,
    ///
    /// Del node probability
    /// `Table[Del, node]`
    ///
    pub d: NodeVec<M>,
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

impl<const M: usize> PHMMTable<M> {
    pub fn new(
        is_dense: bool,
        n_nodes: usize,
        m: Prob,
        i: Prob,
        d: Prob,
        mb: Prob,
        ib: Prob,
        e: Prob,
    ) -> Self {
        PHMMTable {
            m: NodeVec::new(n_nodes, m, is_dense),
            i: NodeVec::new(n_nodes, i, is_dense),
            d: NodeVec::new(n_nodes, d, is_dense),
            mb,
            ib,
            e,
        }
    }
    pub fn zero(is_dense: bool, n_nodes: usize) -> Self {
        PHMMTable {
            m: NodeVec::new(n_nodes, Prob::zero(), is_dense),
            i: NodeVec::new(n_nodes, Prob::zero(), is_dense),
            d: NodeVec::new(n_nodes, Prob::zero(), is_dense),
            mb: Prob::zero(),
            ib: Prob::zero(),
            e: Prob::zero(),
        }
    }
}

///
/// PHMMTables
///
/// Vec of PHMMTable to store tables of emissions
///
#[derive(Debug, Clone)]
pub struct PHMMTables<const M: usize> {
    pub init_table: PHMMTable<M>,
    pub tables: Vec<PHMMTable<M>>,
    pub is_forward: bool,
}
