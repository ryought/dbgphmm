//!
//! PHMMTable V2
//!
//! Updated version of PHMMTable using SparseVec
//!
use super::sample::State;
use crate::prob::Prob;
use arrayvec::ArrayVec;
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
    pub fn n_active_nodes(&self) -> usize {
        self.m.n_elements()
    }
    pub fn is_dense(&self) -> bool {
        self.m.is_dense()
    }
    ///
    ///
    ///
    pub fn filled_nodes(&self) -> Option<ArrayVec<NodeIndex, MAX_ACTIVE_NODES>> {
        if self.is_dense() {
            None
        } else {
            Some(self.to_nodevec().to_top_k_indexes(self.n_active_nodes()))
        }
    }
    ///
    /// Pick up top-scored nodes
    ///
    pub fn top_nodes(&self, n_nodes: usize) -> ArrayVec<NodeIndex, MAX_ACTIVE_NODES> {
        self.to_nodevec().to_top_k_indexes(n_nodes)
    }
    ///
    /// Pick up top-scored nodes by score ratio
    /// for sorted vector v, pick an element while `log(p[i]/p[0])` is greater than max_ratio.
    ///
    pub fn top_nodes_by_score_ratio(
        &self,
        max_ratio: f64,
    ) -> ArrayVec<NodeIndex, MAX_ACTIVE_NODES> {
        let mut ret = ArrayVec::new();
        let v = self.to_nodevec().to_sorted_arrayvec();
        if !v.is_empty() {
            let (_, p0) = v[0];
            for (node, prob) in v {
                if p0.to_log_value() - prob.to_log_value() < max_ratio {
                    ret.push(node);
                }
            }
        }
        ret
    }
    ///
    /// Pick up top-scored nodes
    ///
    pub fn top_nodes_with_prob(&self, n_nodes: usize) -> Vec<(NodeIndex, Prob)> {
        let v = self.to_nodevec();
        v.to_top_k_indexes(n_nodes)
            .into_iter()
            .map(|i| (i, v[i]))
            .collect()
    }
    ///
    /// Pick up top-scored nodes
    ///
    pub fn top_nodes_with_prob_by_score_ratio(&self, max_ratio: f64) -> Vec<(NodeIndex, Prob)> {
        let v = self.to_nodevec();
        self.top_nodes_by_score_ratio(max_ratio)
            .into_iter()
            .map(|i| (i, v[i]))
            .collect()
    }
    /// Diff of two PHMMTables
    ///
    /// sum of |log(pa)-log(pb)|
    ///
    pub fn log_diff(&self, other: &PHMMTable) -> f64 {
        let cmp = |pa: Prob, pb: Prob| pa.log_diff(pb);
        self.m.diff_by(&other.m, cmp)
            + self.i.diff_by(&other.i, cmp)
            + self.d.diff_by(&other.d, cmp)
            + self.mb.log_diff(other.mb)
            + self.ib.log_diff(other.ib)
            + self.e.log_diff(other.e)
    }
    /// Diff of two PHMMTables
    ///
    /// sum of |pa-pb|
    ///
    pub fn diff(&self, other: &PHMMTable) -> f64 {
        let cmp = |pa: Prob, pb: Prob| pa.diff(pb);
        self.m.diff_by(&other.m, cmp)
            + self.i.diff_by(&other.i, cmp)
            + self.d.diff_by(&other.d, cmp)
            + self.mb.diff(other.mb)
            + self.ib.diff(other.ib)
            + self.e.diff(other.e)
    }
    /// PHMMTable (= three NodeVecs of m/i/d) -> merged NodeVec of m+i+d
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
        self.to_states_normalize(None, false)
    }
    /// with normalized probability
    ///
    ///
    pub fn to_states_normalize(&self, p0: Option<Prob>, is_normalize: bool) -> Vec<(State, Prob)> {
        let z = if is_normalize {
            self.m.sum() + self.i.sum() + self.d.sum()
        } else {
            Prob::one()
        };
        let mut ret = Vec::new();
        for (node, p) in self.m.iter() {
            ret.push((State::Match(node), p / z));
        }
        for (node, p) in self.i.iter() {
            ret.push((State::Ins(node), p / z));
        }
        for (node, p) in self.d.iter() {
            ret.push((State::Del(node), p / z));
        }
        ret.sort_by_key(|(_, p)| *p);
        ret.reverse();

        if let Some(p0) = p0 {
            let mut r = Vec::new();
            let mut p_cum = Prob::zero();
            for (state, p) in ret {
                p_cum += p;
                r.push((state, p));
                if p_cum > p0 {
                    break;
                }
            }
            r
        } else {
            ret
        }
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
                    format!("{}:{}{:.5}", node_info(node), state, prob.to_log_value())
                } else {
                    format!("{}{:.5}", state, prob.to_log_value())
                }
            })
            .join(",")
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
    pub kind: PHMMKind,
}

#[derive(Clone, Debug, Copy, PartialEq)]
pub enum PHMMKind {
    Forward,
    Backward,
    State,
}

impl PHMMTables {
    pub fn n_emissions(&self) -> usize {
        self.tables.len()
    }
    pub fn table(&self, i: usize) -> &PHMMTable {
        &self.tables[i]
    }
    pub fn first_table(&self) -> &PHMMTable {
        self.tables
            .first()
            .expect("cannot get first table because PHMMTables is empty")
    }
    pub fn last_table(&self) -> &PHMMTable {
        self.tables
            .last()
            .expect("cannot get last table because PHMMTables is empty")
    }
    pub fn full_prob(&self) -> Prob {
        match self.kind {
            PHMMKind::Forward => self.last_table().e,
            PHMMKind::Backward => self.first_table().mb,
            _ => panic!(),
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
        // TODO
        // assert!(merged_index <= self.n_emissions() + 1);
        match self.kind {
            PHMMKind::Forward => {
                if merged_index == 0 {
                    &self.init_table
                } else {
                    self.table(merged_index - 1)
                }
            }
            PHMMKind::Backward => {
                if merged_index >= self.n_emissions() {
                    &self.init_table
                } else {
                    self.table(merged_index)
                }
            }
            _ => panic!(),
        }
    }
}

/// Store output of forward and backward of a emission sequence.
///
/// # Struct
///
/// * forward: PHMMTables
/// * backward: PHMMTables
///
///
/// # Methods
///
/// calculations that requires both forward and backward result should be in this section.
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
        assert_eq!(forward.kind, PHMMKind::Forward);
        assert_eq!(backward.kind, PHMMKind::Backward);
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
    ///
    /// Calculate emit probs (the probability of nodes to emit the emission x[i]).
    ///
    /// index uses merged-index, so state prob of emission[0] can be calculated by merged_index=1.
    ///
    pub fn to_emit_probs(&self, merged_index: usize) -> PHMMTable {
        let p = self.to_full_prob_forward();
        let f = self.forward.table_merged(merged_index);
        let b = self.backward.table_merged(merged_index);
        (f * b) / p
    }
    ///
    /// S[i] = F[i] B[i]
    ///
    pub fn to_state_probs_v2(&self) -> PHMMTables {
        PHMMTables {
            init_table: self.to_emit_probs(0),
            tables: (1..=self.n_emissions())
                .map(|i| self.to_emit_probs(i))
                .collect(),
            kind: PHMMKind::State,
        }
    }
}
