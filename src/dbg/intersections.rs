use super::dbg::{Dbg, DbgEdge, DbgNode};
use super::edge_centric::{SimpleEDbg, SimpleEDbgEdge, SimpleEDbgNode};
use crate::kmer::kmer::KmerLike;

///
/// Iterator on de bruijn graph intersections corresponding k-1-mer overlaps
///
/// ## Example (k = 4)
///
/// ATTC ------> TTCG
///
pub struct Intersections<K: KmerLike> {
    edbg: SimpleEDbg<K>,
}

///
/// Intersections
///
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    pub fn intersections(&self) -> Intersections<N::Kmer> {
        let edbg = self.to_edbg();
        unimplemented!();
    }
}
