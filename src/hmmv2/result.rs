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
    ///
    /// This PHMM result represents forward or backward?
    ///
    fn is_forward(&self) -> bool;
    /// The number of emissions that this result stores.
    fn n_emissions(&self) -> usize;
    /// index-access to table
    fn table(&self, index: usize) -> PHMMTableRef;
    /// get init_table
    fn init_table(&self) -> PHMMTableRef;
    fn first_table(&self) -> PHMMTableRef {
        self.table(0)
    }
    fn last_table(&self) -> PHMMTableRef {
        self.table(self.n_emissions() - 1)
    }
    ///
    /// Access to a table by merged_index (0 <= i <= n).
    ///
    /// Forward[i]
    ///  = init_table  if i==0
    ///    table(i-1)  otherwise
    /// Backward[i]
    ///  = init_table  if i==n
    ///    table(i)    otherwise
    ///
    fn table_merged(&self, merged_index: usize) -> PHMMTableRef {
        if self.is_forward() {
            if merged_index == 0 {
                self.init_table()
            } else {
                self.table(merged_index - 1)
            }
        } else {
            if merged_index >= self.n_emissions() {
                self.init_table()
            } else {
                self.table(merged_index)
            }
        }
    }
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
    pub is_forward: bool,
}

impl PHMMResultLike for PHMMResult {
    fn is_forward(&self) -> bool {
        self.is_forward
    }
    fn n_emissions(&self) -> usize {
        self.tables.len()
    }
    fn init_table(&self) -> PHMMTableRef {
        PHMMTableRef::Dense(&self.init_table)
    }
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

impl PHMMResultLike for PHMMResultSparse {
    fn is_forward(&self) -> bool {
        self.is_forward
    }
    fn n_emissions(&self) -> usize {
        self.tables_warmup.len() + self.tables_sparse.len()
    }
    /// TODO getter of table
    /// TODO index should be enum?
    fn table(&self, index: usize) -> PHMMTableRef {
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
    fn init_table(&self) -> PHMMTableRef {
        PHMMTableRef::Dense(&self.init_table)
    }
}
