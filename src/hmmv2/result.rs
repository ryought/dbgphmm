//!
//! PHMMResult definitions
//!
//! * For dense
//! * For sparse
//!
use super::table::PHMMTable;
use crate::prob::Prob;
use crate::vector::{DenseStorage, SparseStorage, Storage};

/// Struct that stores Forward/Backward algorithm result
/// for the given emissions
///
/// the length of the PHMMResult.tables will be
/// equal to the length of emissions
#[derive(Debug, Clone)]
pub struct PHMMResult<S: Storage<Item = Prob>> {
    pub init_table: PHMMTable<S>,
    pub tables: Vec<PHMMTable<S>>,
}

impl<S: Storage<Item = Prob>> PHMMResult<S> {
    /// The number of emissions that this result stores.
    pub fn n_emissions(&self) -> usize {
        self.tables.len()
    }
}

///
/// PHMMResult for sparse calculation
///
pub struct PHMMResultSparse {
    pub init_table: PHMMTable<DenseStorage<Prob>>,
    pub tables_warmup: Vec<PHMMTable<DenseStorage<Prob>>>,
    pub tables_sparse: Vec<PHMMTable<SparseStorage<Prob>>>,
}

impl PHMMResultSparse {
    /// The number of emissions that this result stores.
    pub fn n_emissions(&self) -> usize {
        self.tables_warmup.len() + self.tables_sparse.len()
    }
}
