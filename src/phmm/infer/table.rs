//!
//! PHMMTable
//!

use super::super::State;
use crate::prob::{p, Prob};
use sparsevec::SparseVec;

use itertools::Itertools;
use petgraph::graph::NodeIndex;
// use std::ops::{Add, AddAssign, Div, Mul, MulAssign};

///
/// PHMMTable
///
/// A Vector `T[Node, Type]` to compute forward/backward DP.
/// * Node is normal graph node and End node
/// * Type is Match/Ins/Del
///
/// * m[node], i[node], d[node]: SparseVec<Node, Prob>
/// * mb, ib, e: Prob
///
#[derive(Debug, Clone)]
pub struct PHMMTable {
    /// `Table[Match, node]`
    pub m: NodeVec,
    /// `Table[Ins, node]`
    pub i: NodeVec,
    /// `Table[Del, node]`
    pub d: NodeVec,
    /// `Table[Match, Begin]`
    pub mb: Prob,
    /// `Table[Ins, Begin]`
    pub ib: Prob,
    /// `Table[End]`
    pub e: Prob,
}

///
/// Maximum number of continuous deletion (D -> D transition)
/// allowed in pHMM.
///
pub const MAX_DEL: usize = 4;

///
/// Maximum number of active nodes (= nodes which have high probability)
///
pub const MAX_ACTIVE_NODES: usize = 40;

///
/// SparseVec used in PHMMTable
///
/// `SparseVec<Prob, NodeIndex, MAX_ACTIVE_NODES>`
///
pub type NodeVec = SparseVec<Prob, NodeIndex, MAX_ACTIVE_NODES>;

//
//
// implementations
//
//

///
/// Constructors of PHMMTable
///
impl PHMMTable {
    ///
    /// Construct PHMMTable with Dense vector with specified default values
    ///
    pub fn new_dense(
        n_nodes: usize,
        m: Prob,
        i: Prob,
        d: Prob,
        mb: Prob,
        ib: Prob,
        e: Prob,
    ) -> Self {
        PHMMTable {
            m: NodeVec::new_dense(n_nodes, m),
            i: NodeVec::new_dense(n_nodes, i),
            d: NodeVec::new_dense(n_nodes, d),
            mb,
            ib,
            e,
        }
    }
    ///
    /// Construct PHMMTable with Sparse vector with specified default values
    ///
    pub fn new_sparse(
        n_nodes: usize,
        m: Prob,
        i: Prob,
        d: Prob,
        mb: Prob,
        ib: Prob,
        e: Prob,
    ) -> Self {
        PHMMTable {
            m: NodeVec::new_sparse(n_nodes, m),
            i: NodeVec::new_sparse(n_nodes, i),
            d: NodeVec::new_sparse(n_nodes, d),
            mb,
            ib,
            e,
        }
    }
    ///
    /// Construct PHMMTable with Dense vector all filled with 0
    ///
    pub fn zero_dense(n_nodes: usize) -> Self {
        PHMMTable::new_dense(n_nodes, p(0.0), p(0.0), p(0.0), p(0.0), p(0.0), p(0.0))
    }
    ///
    /// Construct PHMMTable with Sparse vector all filled with 0
    ///
    pub fn zero_sparse(n_nodes: usize) -> Self {
        PHMMTable::new_sparse(n_nodes, p(0.0), p(0.0), p(0.0), p(0.0), p(0.0), p(0.0))
    }
    ///
    /// Number of nodes in graph
    ///
    pub fn n_nodes(&self) -> usize {
        self.m.len()
    }
    ///
    /// NodeVec is dense/sparse?
    ///
    pub fn is_dense(&self) -> bool {
        self.m.is_dense()
    }
    ///
    /// PHMMTable into Vec<(State, Prob)>
    ///
    /// Using in creating summary_string in `PHMMTable::to_summary_string`
    ///
    pub fn to_states(&self) -> Vec<(State, Prob)> {
        let mut ret = Vec::new();
        for (node, p) in self.m.iter() {
            ret.push((State::Match(node), p));
        }
        for (node, p) in self.i.iter() {
            ret.push((State::Ins(node), p));
        }
        for (node, p) in self.d.iter() {
            ret.push((State::Del(node), p));
        }
        ret.sort_by_key(|(_, p)| *p);
        ret.reverse();
        ret
    }
}

/*
impl PHMMTable {
    /// Convert to the nodevec containing the sum of hidden states for each node
    /// `v[node] = t.m[node] + t.i[node] + t.d[node]`
    ///
    pub fn to_nodevec(&self) -> NodeVec<S> {
        let mut v = NodeVec::new(self.n_nodes(), Prob::from_prob(0.0));
        for (node, p_m) in self.m.iter() {
            v[node] += p_m;
        }
        for (node, p_i) in self.i.iter() {
            v[node] += p_i;
        }
        for (node, p_d) in self.d.iter() {
            v[node] += p_d;
        }
        v
    }
}

/// Active nodes related methods
impl<S: Storage<Item = Prob>> PHMMTable<S> {
    /// refresh self.active_nodes by its own probabilities
    pub fn refresh_active_nodes(&mut self, n_active_nodes: usize) {
        self.active_nodes = self.active_nodes_from_prob(n_active_nodes);
    }
    /// Determine latest active_nodes from the current probabilities in the table
    ///
    pub fn active_nodes_from_prob(&self, n_active_nodes: usize) -> ActiveNodes {
        ActiveNodes::from_nodevec(&self.to_nodevec(), n_active_nodes)
    }
}

//
// Diff related functions
//
impl<S: Storage<Item = Prob>> PHMMTable<S> {
    ///
    /// Measureing difference between two phmm tables (M/I/D)
    ///
    /// If differences in MB/IB/E matters, please use `diff_all`
    ///
    pub fn diff<T: Storage<Item = Prob>>(&self, other: &PHMMTable<T>) -> f64 {
        self.m.diff(&other.m) + self.i.diff(&other.i) + self.d.diff(&other.d)
    }
    ///
    /// Measureing difference between two phmm tables (M/I/D/MB/MI/E)
    ///
    pub fn diff_all<T: Storage<Item = Prob>>(&self, other: &PHMMTable<T>) -> f64 {
        self.m.diff(&other.m)
            + self.i.diff(&other.i)
            + self.d.diff(&other.d)
            + self.mb.diff(other.mb)
            + self.ib.diff(other.ib)
            + self.e.diff(other.e)
    }
}
*/

impl std::fmt::Display for PHMMTable {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // Header
        if self.is_dense() {
            write!(f, "dense")?;
        } else {
            write!(f, "sparse")?;
        }
        writeln!(f, "\tMatch\tIns\tDel")?;
        // Begin state
        writeln!(f, "Begin\t{}\t{}", self.mb, self.ib)?;
        // Normal states
        for i in 0..self.n_nodes() {
            let v = NodeIndex::new(i);
            writeln!(f, "\t{}\t{}\t{}", self.m[v], self.i[v], self.d[v])?;
        }
        // End state
        writeln!(f, "End\t{}", self.e)
    }
}

//
// summarized to_string for limited area
//

impl PHMMTable {
    ///
    /// one-line summary with top-10 nodes
    ///
    pub fn to_summary_string<F: Fn(NodeIndex) -> String>(&self, node_info: F) -> String {
        self.to_summary_string_n(10, node_info)
    }
    ///
    /// one-line summary
    ///
    pub fn to_summary_string_n<F: Fn(NodeIndex) -> String>(
        &self,
        n: usize,
        node_info: F,
    ) -> String {
        self.to_states()
            .into_iter()
            .take(n)
            .map(|(state, prob)| {
                if let Some(node) = state.to_node_index() {
                    format!("{}:{}{:.5}", node_info(node), state, prob.to_value())
                } else {
                    format!("{}{:.5}", state, prob.to_value())
                }
            })
            .join(",")
    }
}
