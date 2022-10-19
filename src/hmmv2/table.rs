//!
//! Table definitions
//!
//! ## PHMMTable
//!
//! the prob assigned for each (nodes, node_type)
//!
//! F[Match,v] or B[Match,v]
//!
use super::sample::State;
use crate::graph::active_nodes::ActiveNodes;
use crate::prob::{p, Prob};
use crate::vector::{DenseStorage, NodeVec, SparseStorage, Storage};
pub use petgraph::graph::NodeIndex;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign};

/// Maximum number of continuous deletion (D -> D transition)
/// allowed in pHMM.
pub const MAX_DEL: usize = 4;

/// Struct for storing Forward/Backward intermediate result
/// for an emission.
///
/// Corresponds to a vector `T[node, type]`
/// `node` is either normal or begin or end node.
/// `type` is either Match, Ins, Del.
#[derive(Debug, Clone, PartialEq)]
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
    ///
    pub active_nodes: ActiveNodes,
}

//
// PHMMTables
//

/// Constructors of PHMMTable
impl<S: Storage<Item = Prob>> PHMMTable<S> {
    pub fn new(n_nodes: usize, m: Prob, i: Prob, d: Prob, mb: Prob, ib: Prob, e: Prob) -> Self {
        PHMMTable {
            m: NodeVec::new(n_nodes, m),
            i: NodeVec::new(n_nodes, i),
            d: NodeVec::new(n_nodes, d),
            mb,
            ib,
            e,
            active_nodes: ActiveNodes::All,
        }
    }
    pub fn new_with_active_nodes(
        n_nodes: usize,
        m: Prob,
        i: Prob,
        d: Prob,
        mb: Prob,
        ib: Prob,
        e: Prob,
        active_nodes: ActiveNodes,
    ) -> Self {
        PHMMTable {
            m: NodeVec::new(n_nodes, m),
            i: NodeVec::new(n_nodes, i),
            d: NodeVec::new(n_nodes, d),
            mb,
            ib,
            e,
            active_nodes,
        }
    }
    pub fn zero(n_nodes: usize) -> Self {
        PHMMTable {
            m: NodeVec::new(n_nodes, Prob::from_prob(0.0)),
            i: NodeVec::new(n_nodes, Prob::from_prob(0.0)),
            d: NodeVec::new(n_nodes, Prob::from_prob(0.0)),
            mb: Prob::from_prob(0.0),
            ib: Prob::from_prob(0.0),
            e: Prob::from_prob(0.0),
            active_nodes: ActiveNodes::All,
        }
    }
    pub fn zero_with_active_nodes(n_nodes: usize, active_nodes: ActiveNodes) -> Self {
        PHMMTable {
            m: NodeVec::new(n_nodes, Prob::from_prob(0.0)),
            i: NodeVec::new(n_nodes, Prob::from_prob(0.0)),
            d: NodeVec::new(n_nodes, Prob::from_prob(0.0)),
            mb: Prob::from_prob(0.0),
            ib: Prob::from_prob(0.0),
            e: Prob::from_prob(0.0),
            active_nodes,
        }
    }
}

/// Accessors of PHMMTable
impl<S: Storage<Item = Prob>> PHMMTable<S> {
    /// Get the number of nodes in the PHMM Table
    pub fn n_nodes(&self) -> usize {
        // self.i and self.d has also the same length
        self.m.len()
    }
    /// Convert to the DenseStorage-backend PHMMTable
    pub fn to_dense(&self) -> PHMMTable<DenseStorage<Prob>> {
        PHMMTable {
            m: self.m.to_dense(),
            i: self.i.to_dense(),
            d: self.d.to_dense(),
            mb: self.mb,
            ib: self.ib,
            e: self.e,
            active_nodes: self.active_nodes.clone(),
        }
    }
    /// Convert to the SparseStorage-backend PHMMTable
    pub fn to_sparse(&self, default_value: Prob) -> PHMMTable<SparseStorage<Prob>> {
        PHMMTable {
            m: self.m.to_sparse(default_value),
            i: self.i.to_sparse(default_value),
            d: self.d.to_sparse(default_value),
            mb: self.mb,
            ib: self.ib,
            e: self.e,
            active_nodes: self.active_nodes.clone(),
        }
    }
    /// Convert to the SparseStorage-backend PHMMTable
    /// by taking only active_nodes
    pub fn to_sparse_active_nodes(&self, n_active_nodes: usize) -> PHMMTable<SparseStorage<Prob>> {
        let active_nodes = self.active_nodes_from_prob(n_active_nodes);
        match active_nodes {
            ActiveNodes::Only(nodes) => PHMMTable {
                m: self.m.to_sparse_by_indexes(p(0.0), &nodes),
                i: self.i.to_sparse_by_indexes(p(0.0), &nodes),
                d: self.d.to_sparse_by_indexes(p(0.0), &nodes),
                mb: self.mb,
                ib: self.ib,
                e: self.e,
                active_nodes: ActiveNodes::Only(nodes),
            },
            ActiveNodes::All => unreachable!(),
        }
    }
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
    /// Measureing difference between two phmm tables
    ///
    pub fn diff<T: Storage<Item = Prob>>(&self, other: &PHMMTable<T>) -> f64 {
        self.m.diff(&other.m) + self.i.diff(&other.i) + self.d.diff(&other.d)
    }
}

/*
/// for approx `assert_abs_diff_eq`
use approx::AbsDiffEq;
impl<S: Storage<Item = Prob>> AbsDiffEq for PHMMTable<S> {
    type Epsilon = <<S as Storage>::Item as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        S::Item::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let m = NodeVec::abs_diff_eq(&self.m, &other.m, epsilon);
        let i = NodeVec::abs_diff_eq(&self.i, &other.i, epsilon);
        let d = NodeVec::abs_diff_eq(&self.d, &other.d, epsilon);
        let ib = Prob::abs_diff_eq(&self.ib, &other.ib, epsilon);
        let mb = Prob::abs_diff_eq(&self.mb, &other.mb, epsilon);
        println!("{} {} {} {} {}", m, i, d, ib, mb);
        m && i && d
    }
}
*/

impl<S: Storage<Item = Prob>> std::fmt::Display for PHMMTable<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // Header
        if S::is_dense() {
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
            match &self.active_nodes {
                ActiveNodes::Only(nodes) => {
                    if nodes.contains(&v) {
                        write!(f, "{}+", i)?;
                    } else {
                        write!(f, "{}", i)?;
                    }
                }
                ActiveNodes::All => {
                    write!(f, "{}*", i)?;
                }
            };
            writeln!(f, "\t{}\t{}\t{}", self.m[v], self.i[v], self.d[v])?;
        }
        // End state
        writeln!(f, "End\t{}", self.e)
    }
}

// Add
impl<'a, 'b, Sa, Sb> Add<&'a PHMMTable<Sa>> for &'b PHMMTable<Sb>
where
    Sa: Storage<Item = Prob>,
    Sb: Storage<Item = Prob>,
{
    type Output = PHMMTable<Sb>;
    fn add(self, other: &'a PHMMTable<Sa>) -> Self::Output {
        PHMMTable {
            m: &self.m + &other.m,
            i: &self.i + &other.i,
            d: &self.d + &other.d,
            mb: self.mb + other.mb,
            ib: self.ib + other.ib,
            e: self.e + other.e,
            active_nodes: ActiveNodes::All,
        }
    }
}

// AddAssign
impl<'a, S> AddAssign<&'a PHMMTable<S>> for PHMMTable<S>
where
    S: Storage<Item = Prob>,
{
    fn add_assign(&mut self, other: &'a PHMMTable<S>) {
        self.m += &other.m;
        self.i += &other.i;
        self.d += &other.d;
        self.mb += other.mb;
        self.ib += other.ib;
        self.e += other.e;
    }
}

// Mul
impl<'a, 'b, Sa, Sb> Mul<&'a PHMMTable<Sa>> for &'b PHMMTable<Sb>
where
    Sa: Storage<Item = Prob>,
    Sb: Storage<Item = Prob>,
{
    type Output = PHMMTable<Sb>;
    fn mul(self, other: &'a PHMMTable<Sa>) -> Self::Output {
        PHMMTable {
            m: &self.m * &other.m,
            i: &self.i * &other.i,
            d: &self.d * &other.d,
            mb: self.mb * other.mb,
            ib: self.ib * other.ib,
            e: self.e * other.e,
            active_nodes: ActiveNodes::All,
        }
    }
}

// MulAssign
impl<'a, S> MulAssign<&'a PHMMTable<S>> for PHMMTable<S>
where
    S: Storage<Item = Prob>,
{
    fn mul_assign(&mut self, other: &'a PHMMTable<S>) {
        self.m *= &other.m;
        self.i *= &other.i;
        self.d *= &other.d;
        self.mb *= other.mb;
        self.ib *= other.ib;
        self.e *= other.e;
    }
}

// Div by constant
impl<S> Div<Prob> for PHMMTable<S>
where
    S: Storage<Item = Prob>,
{
    type Output = PHMMTable<S>;
    fn div(self, other: Prob) -> Self::Output {
        PHMMTable {
            m: self.m / other,
            i: self.i / other,
            d: self.d / other,
            mb: self.mb / other,
            ib: self.ib / other,
            e: self.e / other,
            active_nodes: ActiveNodes::All,
        }
    }
}

// Sum
impl<S> std::iter::Sum for PHMMTable<S>
where
    S: Storage<Item = Prob>,
{
    fn sum<I>(iter: I) -> PHMMTable<S>
    where
        I: Iterator<Item = PHMMTable<S>>,
    {
        iter.reduce(|mut a, b| {
            a += &b;
            a
        })
        .unwrap()
    }
}

// Product
impl<S> std::iter::Product for PHMMTable<S>
where
    S: Storage<Item = Prob>,
{
    fn product<I>(iter: I) -> PHMMTable<S>
    where
        I: Iterator<Item = PHMMTable<S>>,
    {
        iter.reduce(|mut a, b| {
            a *= &b;
            a
        })
        .unwrap()
    }
}

//
// PHMMResults
//

// TODO add display for PHMMResult

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;
    use crate::prob::p;
    use crate::vector::DenseStorage;

    ///
    /// Helper function for PHMMTable tests
    /// check if `|t1[i] - t2[i]|` is small for all `i`.
    ///
    fn check_is_similar_table<S: Storage<Item = Prob>>(t1: &PHMMTable<S>, t2: &PHMMTable<S>) {
        assert_eq!(t1.n_nodes(), t2.n_nodes());
        let max_diff = 0.0000001;
        // check t.m, t.i, t.d
        for i in 0..t1.n_nodes() {
            let v = NodeIndex::new(i);
            assert!((t1.m[v].to_value() - t2.m[v].to_value()).abs() < max_diff);
            assert!((t1.i[v].to_value() - t2.i[v].to_value()).abs() < max_diff);
            assert!((t1.d[v].to_value() - t2.d[v].to_value()).abs() < max_diff);
        }
        // check t.mb, t.ib, t.e
        assert!((t1.mb.to_value() - t2.mb.to_value()).abs() < max_diff);
        assert!((t1.ib.to_value() - t2.ib.to_value()).abs() < max_diff);
        assert!((t1.e.to_value() - t2.e.to_value()).abs() < max_diff);
    }

    #[test]
    fn hmm_table_add() {
        let mut t1: PHMMTable<DenseStorage<Prob>> = PHMMTable::zero(5);
        let mut t2: PHMMTable<DenseStorage<Prob>> = PHMMTable::zero(5);
        t1.m[NodeIndex::new(0)] = Prob::from_prob(0.5);
        t2.m[NodeIndex::new(0)] = Prob::from_prob(0.5);
        t1.d[NodeIndex::new(0)] = Prob::from_prob(0.3);
        t2.d[NodeIndex::new(0)] = Prob::from_prob(0.3);
        t1.i[NodeIndex::new(1)] = Prob::from_prob(0.2);
        t2.i[NodeIndex::new(1)] = Prob::from_prob(0.0);
        t1.e = Prob::from_prob(0.2);
        t2.e = Prob::from_prob(0.4);
        t1.mb = Prob::from_prob(0.01);
        t2.mb = Prob::from_prob(0.02);
        println!("{}", t1);
        println!("{}", t2);

        let mut t_add: PHMMTable<DenseStorage<Prob>> = PHMMTable::zero(5);
        t_add.m[NodeIndex::new(0)] = Prob::from_prob(1.0);
        t_add.d[NodeIndex::new(0)] = Prob::from_prob(0.6);
        t_add.i[NodeIndex::new(1)] = Prob::from_prob(0.2);
        t_add.e = Prob::from_prob(0.6);
        t_add.mb = Prob::from_prob(0.03);
        println!("t_add\n{}", t_add);
        let t_add2 = &t1 + &t2;
        println!("{}", t_add2);
        check_is_similar_table(&t_add, &t_add2);

        // add assign
        let mut t3 = t1.clone();
        t3 += &t2;
        check_is_similar_table(&t_add, &t3);

        // sum
        let ts: Vec<PHMMTable<DenseStorage<Prob>>> = vec![t1, t2];
        let sum: PHMMTable<DenseStorage<Prob>> = ts.into_iter().sum();
        println!("sum\n{}", sum);
        check_is_similar_table(&t_add, &sum);
    }
    #[test]
    fn hmm_table_mul() {
        let mut t1: PHMMTable<DenseStorage<Prob>> = PHMMTable::zero(5);
        let mut t2: PHMMTable<DenseStorage<Prob>> = PHMMTable::zero(5);
        t1.m[NodeIndex::new(0)] = Prob::from_prob(0.5);
        t2.m[NodeIndex::new(0)] = Prob::from_prob(0.5);
        t1.d[NodeIndex::new(0)] = Prob::from_prob(0.3);
        t2.d[NodeIndex::new(0)] = Prob::from_prob(0.3);
        t1.i[NodeIndex::new(1)] = Prob::from_prob(0.2);
        t2.i[NodeIndex::new(1)] = Prob::from_prob(0.0);
        t1.m[NodeIndex::new(1)] = Prob::from_prob(1.0);
        t2.m[NodeIndex::new(1)] = Prob::from_prob(0.5);
        t1.e = Prob::from_prob(0.2);
        t2.e = Prob::from_prob(0.4);
        t1.mb = Prob::from_prob(0.01);
        t2.mb = Prob::from_prob(0.02);
        println!("{}", t1);
        println!("{}", t2);

        let mut t_mul: PHMMTable<DenseStorage<Prob>> = PHMMTable::zero(5);
        t_mul.m[NodeIndex::new(0)] = Prob::from_prob(0.25);
        t_mul.d[NodeIndex::new(0)] = Prob::from_prob(0.09);
        t_mul.i[NodeIndex::new(1)] = Prob::from_prob(0.0);
        t_mul.m[NodeIndex::new(1)] = Prob::from_prob(0.5);
        t_mul.e = Prob::from_prob(0.08);
        t_mul.mb = Prob::from_prob(0.0002);
        let t_mul2 = &t1 * &t2;
        check_is_similar_table(&t_mul, &t_mul2);

        // mul asign
        let mut t3 = t1.clone();
        t3 *= &t2;
        check_is_similar_table(&t_mul, &t3);

        // product
        let ts: Vec<PHMMTable<DenseStorage<Prob>>> = vec![t1, t2];
        let prd: PHMMTable<DenseStorage<Prob>> = ts.into_iter().product();
        println!("t_mul\n{}", t_mul);
        println!("prd\n{}", prd);
        // TODO
        check_is_similar_table(&t_mul, &prd);
    }
    #[test]
    fn hmm_table_to_nodevec() {
        let mut t1: PHMMTable<DenseStorage<Prob>> = PHMMTable::zero(5);
        t1.m[ni(0)] = p(0.5);
        t1.d[ni(0)] = p(0.01);
        t1.i[ni(0)] = p(0.02);
        t1.i[ni(1)] = p(0.2);
        t1.m[ni(1)] = p(1.0);
        t1.e = p(0.2);
        t1.mb = p(0.01);
        let v = t1.to_nodevec();
        println!("{:?}", v);
        assert_abs_diff_eq!(v[ni(0)], p(0.5) + p(0.01) + p(0.02));
        assert_abs_diff_eq!(v[ni(1)], p(0.2) + p(1.0));
        assert!(v[ni(2)].is_zero());
        assert!(v[ni(3)].is_zero());

        let mut t1: PHMMTable<SparseStorage<Prob>> = PHMMTable::zero(5);
        t1.m[ni(0)] = p(0.5);
        t1.d[ni(0)] = p(0.01);
        t1.i[ni(0)] = p(0.02);
        t1.i[ni(1)] = p(0.2);
        t1.m[ni(1)] = p(1.0);
        t1.e = p(0.2);
        t1.mb = p(0.01);
        let v = t1.to_nodevec();
        println!("{:?}", v);
        assert_abs_diff_eq!(v[ni(0)], p(0.5) + p(0.01) + p(0.02));
        assert_abs_diff_eq!(v[ni(1)], p(0.2) + p(1.0));
        assert!(v[ni(2)].is_zero());
        assert!(v[ni(3)].is_zero());
    }
    #[test]
    fn hmm_table_refresh() {
        let mut t: PHMMTable<DenseStorage<Prob>> = PHMMTable::zero(5);
        t.m[ni(0)] = p(0.5);
        t.d[ni(0)] = p(0.01);
        t.i[ni(0)] = p(0.02);
        t.i[ni(1)] = p(0.2);
        t.m[ni(1)] = p(1.0);
        t.refresh_active_nodes(1);
        println!("{}", t);
        assert_eq!(t.active_nodes, ActiveNodes::Only(vec![ni(1)]))
    }
    #[test]
    fn hmm_table_diff() {
        let mut t1: PHMMTable<DenseStorage<Prob>> = PHMMTable::zero(5);
        t1.m[ni(0)] = p(0.5);
        t1.d[ni(0)] = p(0.01);
        t1.i[ni(0)] = p(0.02);
        t1.i[ni(1)] = p(0.2);
        t1.m[ni(1)] = p(1.0);
        let mut t2: PHMMTable<DenseStorage<Prob>> = PHMMTable::zero(5);
        t2.m[ni(0)] = p(0.51);
        t2.d[ni(0)] = p(0.02);
        t2.i[ni(0)] = p(0.02);
        t2.i[ni(1)] = p(0.21);
        t2.m[ni(1)] = p(1.0);
        println!("d(t1,t2)={}", t1.diff(&t2));
        println!("d(t2,t1)={}", t2.diff(&t1));
        assert_abs_diff_eq!(t1.diff(&t2), 0.03);
        assert_abs_diff_eq!(t2.diff(&t1), 0.03);
        // use approx
        println!("{}", t1);
        println!("{}", t2);
        // assert!(abs_diff_eq!(t1, t1));
        // assert!(abs_diff_eq!(t1, t1, epsilon = 1.0));
    }
}
