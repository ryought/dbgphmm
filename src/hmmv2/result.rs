//!
//! PHMMResult definitions
//!
//! * For dense
//! * For sparse
//!
use super::table::PHMMTable;
use crate::prob::Prob;
use crate::vector::Storage;

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
