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
    ///
    pub active_nodes: ActiveNodes,
}

#[derive(Debug, Clone)]
pub enum ActiveNodes {
    /// All nodes are active
    All,
    /// Only nodes specified are active
    Only(Vec<NodeIndex>),
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
}

impl<S: Storage<Item = Prob>> std::fmt::Display for PHMMTable<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // Header
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
        self.mb = self.mb + other.mb;
        self.ib = self.ib + other.ib;
        self.e = self.e + other.e;
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
        self.mb = self.mb * other.mb;
        self.ib = self.ib * other.ib;
        self.e = self.e * other.e;
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
}
