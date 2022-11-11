//!
//! Constructor of draft dbg
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::common::{CopyNum, Seq};

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    pub fn create_draft_from_seqs<T>(k: usize, seqs: T) -> Self
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        unimplemented!();
    }
}
