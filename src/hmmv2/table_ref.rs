//!
//! Table enums
//!
use super::table::PHMMTable;
use crate::prob::Prob;
use crate::vector::{DenseStorage, SparseStorage, Storage};
pub use petgraph::graph::NodeIndex;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign};

///
///
///
pub enum PHMMTableRef<'a> {
    Dense(&'a PHMMTable<DenseStorage<Prob>>),
    Sparse(&'a PHMMTable<SparseStorage<Prob>>),
}

///
///
///
pub enum PHMMTableEnum {
    Dense(PHMMTable<DenseStorage<Prob>>),
    Sparse(PHMMTable<SparseStorage<Prob>>),
}

impl<'a> PHMMTableRef<'a> {
    /// Get a value from `PHMMTable.m` by node index
    pub fn m(&self, node: NodeIndex) -> Prob {
        match self {
            PHMMTableRef::Dense(t) => t.m[node],
            PHMMTableRef::Sparse(t) => t.m[node],
        }
    }
    /// Get a value from `PHMMTable.i` by node index
    pub fn i(&self, node: NodeIndex) -> Prob {
        match self {
            PHMMTableRef::Dense(t) => t.i[node],
            PHMMTableRef::Sparse(t) => t.i[node],
        }
    }
    /// Get a value from `PHMMTable.d` by node index
    pub fn d(&self, node: NodeIndex) -> Prob {
        match self {
            PHMMTableRef::Dense(t) => t.d[node],
            PHMMTableRef::Sparse(t) => t.d[node],
        }
    }
    /// check if backend is dense or sparse?
    pub fn is_dense(&self) -> bool {
        match self {
            PHMMTableRef::Dense(_) => true,
            PHMMTableRef::Sparse(_) => false,
        }
    }
    /// convert to dense backed table
    pub fn to_dense(&self) -> PHMMTable<DenseStorage<Prob>> {
        match self {
            PHMMTableRef::Dense(t) => t.to_dense(),
            PHMMTableRef::Sparse(t) => t.to_dense(),
        }
    }
}

//
// Add
//
impl<'a, 'b> Add<&'a PHMMTableRef<'a>> for &'b PHMMTableRef<'b> {
    type Output = PHMMTable<DenseStorage<Prob>>;
    fn add(self, other: &'a PHMMTableRef) -> Self::Output {
        match (self, other) {
            // if dense is involved, return should be dense
            // and add sparse to dense, not dense to sparse.
            (&PHMMTableRef::Dense(d1), &PHMMTableRef::Dense(d2)) => d1 + d2,
            (&PHMMTableRef::Dense(d), &PHMMTableRef::Sparse(s)) => d + s,
            (&PHMMTableRef::Sparse(s), &PHMMTableRef::Dense(d)) => d + s,
            // TODO if sparse-vec only, returned-vec can be sparse
            // for now, all returned-vec is dense
            (&PHMMTableRef::Sparse(s1), &PHMMTableRef::Sparse(s2)) => &s1.to_dense() + s2,
        }
    }
}

//
// Mul
//
impl<'a, 'b> Mul<&'a PHMMTableRef<'a>> for &'b PHMMTableRef<'b> {
    type Output = PHMMTable<DenseStorage<Prob>>;
    fn mul(self, other: &'a PHMMTableRef) -> Self::Output {
        match (self, other) {
            // if dense is involved, return should be dense
            // and add sparse to dense, not dense to sparse.
            (&PHMMTableRef::Dense(d1), &PHMMTableRef::Dense(d2)) => d1 * d2,
            (&PHMMTableRef::Dense(d), &PHMMTableRef::Sparse(s)) => d * s,
            (&PHMMTableRef::Sparse(s), &PHMMTableRef::Dense(d)) => d * s,
            // TODO if sparse-vec only, returned-vec can be sparse
            // for now, all returned-vec is dense
            (&PHMMTableRef::Sparse(s1), &PHMMTableRef::Sparse(s2)) => &s1.to_dense() * s2,
        }
    }
}

//
// Diff
//
impl<'a> PHMMTableRef<'a> {
    ///
    /// Measureing difference between two phmm tables
    ///
    pub fn diff(&self, other: &PHMMTableRef) -> f64 {
        match (self, other) {
            (&PHMMTableRef::Dense(s), &PHMMTableRef::Dense(o)) => s.diff(o),
            (&PHMMTableRef::Dense(s), &PHMMTableRef::Sparse(o)) => s.diff(o),
            (&PHMMTableRef::Sparse(s), &PHMMTableRef::Dense(o)) => s.diff(o),
            (&PHMMTableRef::Sparse(s), &PHMMTableRef::Sparse(o)) => s.diff(o),
        }
    }
}

//
// Display
//
impl<'a> std::fmt::Display for PHMMTableRef<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            PHMMTableRef::Dense(t) => write!(f, "{}", t),
            PHMMTableRef::Sparse(t) => write!(f, "{}", t),
        }
    }
}

//
// tests
//
#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;
    use crate::prob::p;

    #[test]
    fn hmmv2_table_ref_add_mul() {
        let mut t1: PHMMTable<DenseStorage<Prob>> = PHMMTable::zero(5);
        t1.m[ni(0)] = p(0.5);
        t1.d[ni(0)] = p(0.01);
        t1.i[ni(0)] = p(0.02);
        t1.i[ni(1)] = p(0.2);
        t1.m[ni(1)] = p(1.0);
        let mut t2: PHMMTable<SparseStorage<Prob>> = PHMMTable::zero(5);
        t2.m[ni(0)] = p(0.5);
        t2.d[ni(0)] = p(0.01);
        t2.i[ni(0)] = p(0.02);
        t2.i[ni(1)] = p(0.2);
        t2.m[ni(1)] = p(1.0);
        let r1 = PHMMTableRef::Dense(&t1);
        let r2 = PHMMTableRef::Sparse(&t2);
        assert!(r1.is_dense());
        assert!(!r2.is_dense());
        println!("{}", &t1 * &t2);
        println!("{}", &r1 * &r2);
        // TODO
        // assert_eq!(&t1 * &t2, &r1 * &r2);
    }
}
