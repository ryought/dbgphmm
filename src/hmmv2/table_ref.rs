//!
//! Table enums
//!
use super::table::PHMMTable;
use crate::prob::Prob;
use crate::vector::{DenseStorage, SparseStorage, Storage};
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

//
// Add
//
impl<'a, 'b> Add<&'a PHMMTableRef<'a>> for &'b PHMMTableRef<'b> {
    type Output = PHMMTable<DenseStorage<Prob>>;
    fn add(self, other: &'a PHMMTableRef) -> Self::Output {
        match (self, other) {
            // if dense is involved, return should be dense
            (&PHMMTableRef::Sparse(s), &PHMMTableRef::Dense(o)) => o + s,
            (&PHMMTableRef::Dense(s), &PHMMTableRef::Dense(o)) => s + o,
            (&PHMMTableRef::Dense(s), &PHMMTableRef::Sparse(o)) => s + o,
            // TODO if sparse-vec only, returned-vec can be sparse
            // for now, all returned-vec is dense
            (&PHMMTableRef::Sparse(s), &PHMMTableRef::Sparse(o)) => &s.to_dense() + o,
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
            (&PHMMTableRef::Sparse(s), &PHMMTableRef::Dense(o)) => o * s,
            (&PHMMTableRef::Dense(s), &PHMMTableRef::Dense(o)) => s * o,
            (&PHMMTableRef::Dense(s), &PHMMTableRef::Sparse(o)) => s * o,
            // TODO if sparse-vec only, returned-vec can be sparse
            // for now, all returned-vec is dense
            (&PHMMTableRef::Sparse(s), &PHMMTableRef::Sparse(o)) => &s.to_dense() * o,
        }
    }
}
