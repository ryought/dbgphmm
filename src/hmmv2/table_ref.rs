//!
//! Table enums
//!
use super::table::PHMMTable;
use crate::prob::Prob;
use crate::vector::{DenseStorage, SparseStorage, Storage};

///
///
///
pub enum PHMMTableRef<'a> {
    Dense(&'a PHMMTable<DenseStorage<Prob>>),
    Sparse(&'a PHMMTable<SparseStorage<Prob>>),
}
