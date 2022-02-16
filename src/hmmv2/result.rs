//!
//! PHMMResult definitions
//!
//! * For dense
//! * For sparse
//!
use super::table::PHMMTable;
use super::table_ref::PHMMTableRef;
use crate::prob::Prob;
use crate::vector::{DenseStorage, SparseStorage, Storage};

///
/// Trait that generalizes PHMMResultFull and PHMMResultSparse
///
pub trait PHMMResultLike {
    fn n_emissions(&self) -> usize;
    fn table(&self, index: usize) -> PHMMTableRef;
    fn init_table(&self) -> PHMMTableRef;
}

/// Struct that stores Forward/Backward algorithm result
/// for the given emissions
///
/// the length of the PHMMResult.tables will be
/// equal to the length of emissions
#[derive(Debug, Clone)]
pub struct PHMMResult {
    pub init_table: PHMMTable<DenseStorage<Prob>>,
    pub tables: Vec<PHMMTable<DenseStorage<Prob>>>,
}

impl PHMMResultLike for PHMMResult {
    /// The number of emissions that this result stores.
    fn n_emissions(&self) -> usize {
        self.tables.len()
    }
    /// get init_table
    fn init_table(&self) -> PHMMTableRef {
        PHMMTableRef::Dense(&self.init_table)
    }
    /// index-access to table
    fn table(&self, index: usize) -> PHMMTableRef {
        PHMMTableRef::Dense(&self.tables[index])
    }
}

///
/// PHMMResult for sparse calculation
///
pub struct PHMMResultSparse {
    pub init_table: PHMMTable<DenseStorage<Prob>>,
    pub tables_warmup: Vec<PHMMTable<DenseStorage<Prob>>>,
    pub tables_sparse: Vec<PHMMTable<SparseStorage<Prob>>>,
    pub is_forward: bool,
}

impl PHMMResultSparse {
    /// The number of emissions that this result stores.
    pub fn n_emissions(&self) -> usize {
        self.tables_warmup.len() + self.tables_sparse.len()
    }
    /// TODO getter of table
    /// TODO index should be enum?
    pub fn table(&self, index: usize) -> PHMMTableRef {
        if self.is_forward {
            let n_warmup = self.tables_warmup.len();
            if index < n_warmup {
                PHMMTableRef::Dense(&self.tables_warmup[index])
            } else {
                PHMMTableRef::Sparse(&self.tables_sparse[index - n_warmup])
            }
        } else {
            let n_sparse = self.tables_sparse.len();
            if index < n_sparse {
                PHMMTableRef::Sparse(&self.tables_sparse[index])
            } else {
                PHMMTableRef::Dense(&self.tables_warmup[index - n_sparse])
            }
        }
    }
    pub fn init_table(&self) -> PHMMTableRef {
        PHMMTableRef::Dense(&self.init_table)
    }
}
