//!
//! PHMMTable V2
//!
//! Updated version of PHMMTable using SparseVec
//!
use super::super::State;
use crate::prob::{p, Prob};
use itertools::Itertools;
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
pub type NodeVec = SparseVec<Prob, NodeIndex, MAX_ACTIVE_NODES>;

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
pub struct PHMMTable {
    ///
    /// Match node probability
    /// `Table[Match, node]`
    ///
    pub m: NodeVec,
    ///
    /// Ins node probability
    /// `Table[Ins, node]`
    ///
    pub i: NodeVec,
    ///
    /// Del node probability
    /// `Table[Del, node]`
    ///
    pub d: NodeVec,
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

impl PHMMTable {
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
    pub fn n_nodes(&self) -> usize {
        self.m.len()
    }
    pub fn is_dense(&self) -> bool {
        self.m.is_dense()
    }
    /// Pick up top-scored nodes
    ///
    pub fn top_nodes(&self, n_nodes: usize) -> Vec<NodeIndex> {
        unimplemented!();
    }
    pub fn diff(&self, other: &PHMMTable) -> f64 {
        unimplemented!();
        // self.m.diff(&other.m)
        //     + self.i.diff(&other.i)
        //     + self.d.diff(&other.d)
        //     + self.mb.diff(other.mb)
        //     + self.ib.diff(other.ib)
        //     + self.e.diff(other.e)
    }
    ///
    ///
    ///
    pub fn to_nodevec(&self) -> NodeVec {
        let mut v = NodeVec::new(self.n_nodes(), Prob::from_prob(0.0), self.is_dense());
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
    ///
    ///
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
            writeln!(f, "{}\t{}\t{}\t{}", i, self.m[v], self.i[v], self.d[v])?;
        }
        // End state
        writeln!(f, "End\t{}", self.e)
    }
}

impl<'a> std::ops::AddAssign<&'a PHMMTable> for PHMMTable {
    fn add_assign(&mut self, other: &'a PHMMTable) {
        self.m += &other.m;
        self.i += &other.i;
        self.d += &other.d;
        self.mb += other.mb;
        self.ib += other.ib;
        self.e += other.e;
    }
}

impl<'a, 'b> std::ops::Mul<&'a PHMMTable> for &'b PHMMTable {
    type Output = PHMMTable;
    fn mul(self, other: &'a PHMMTable) -> Self::Output {
        PHMMTable {
            m: &self.m * &other.m,
            i: &self.i * &other.i,
            d: &self.d * &other.d,
            mb: self.mb * other.mb,
            ib: self.ib * other.ib,
            e: self.e * other.e,
        }
    }
}

impl std::ops::Div<Prob> for PHMMTable {
    type Output = PHMMTable;
    fn div(self, other: Prob) -> Self::Output {
        PHMMTable {
            m: self.m / other,
            i: self.i / other,
            d: self.d / other,
            mb: self.mb / other,
            ib: self.ib / other,
            e: self.e / other,
        }
    }
}

impl std::iter::Sum for PHMMTable {
    fn sum<I>(iter: I) -> PHMMTable
    where
        I: Iterator<Item = PHMMTable>,
    {
        iter.reduce(|mut a, b| {
            a += &b;
            a
        })
        .unwrap()
    }
}

///
/// PHMMTables
///
/// Vec of PHMMTable to store tables of emissions
///
#[derive(Debug, Clone)]
pub struct PHMMTables {
    pub init_table: PHMMTable,
    pub tables: Vec<PHMMTable>,
    pub is_forward: bool,
}

impl PHMMTables {
    pub fn n_emissions(&self) -> usize {
        self.tables.len()
    }
    pub fn table(&self, i: usize) -> &PHMMTable {
        &self.tables[i]
    }
    pub fn first_table(&self) -> &PHMMTable {
        self.tables.first().unwrap()
    }
    pub fn last_table(&self) -> &PHMMTable {
        self.tables.last().unwrap()
    }
    pub fn is_forward(&self) -> bool {
        self.is_forward
    }
    pub fn full_prob(&self) -> Prob {
        if self.is_forward() {
            self.last_table().e
        } else {
            self.first_table().mb
        }
    }
    ///
    /// Access to a table by merged_index (0 <= i <= n).
    ///
    /// Forward[i]
    ///  = init_table  if i==0
    ///    table(i-1)  otherwise
    ///  = P(emits x[:i] and now in a state)
    /// Backward[i]
    ///  = init_table  if i==n
    ///    table(i)    otherwise
    ///  = P(emits x[i:] | starts from a state)
    ///
    pub fn table_merged(&self, merged_index: usize) -> &PHMMTable {
        assert!(merged_index <= self.n_emissions());
        if self.is_forward() {
            if merged_index == 0 {
                &self.init_table
            } else {
                self.table(merged_index - 1)
            }
        } else {
            if merged_index >= self.n_emissions() {
                &self.init_table
            } else {
                self.table(merged_index)
            }
        }
    }
}

///
///
///
///
#[derive(Debug, Clone)]
pub struct PHMMOutput {
    ///
    /// forward
    ///
    pub forward: PHMMTables,
    ///
    /// backward
    ///
    pub backward: PHMMTables,
}

///
/// Constructors
///
impl PHMMOutput {
    pub fn new(forward: PHMMTables, backward: PHMMTables) -> Self {
        // check forward/backward is created by phmm.forward/backward()
        assert!(forward.is_forward());
        assert!(!backward.is_forward());
        assert_eq!(forward.n_emissions(), backward.n_emissions());
        PHMMOutput { forward, backward }
    }
    pub fn n_emissions(&self) -> usize {
        self.forward.n_emissions()
    }
    /// Calculate the full probability `P(x)` of the given emission `x`
    /// from **forward** result.
    ///
    /// ```text
    /// P(x) = fe_n-1 = P(emits x[0],...,x[n-1] and now in `e` (end state))
    /// ```
    ///
    pub fn to_full_prob_forward(&self) -> Prob {
        self.forward.full_prob()
    }
    /// Calculate the full probability `P(x)` of the given emission `x`
    /// from **backward** result.
    ///
    /// ```text
    /// P(x) = bm_0[b] = P(emits x[0:] | starts from m_b)
    /// ```
    ///
    pub fn to_full_prob_backward(&self) -> Prob {
        self.backward.full_prob()
    }
}
