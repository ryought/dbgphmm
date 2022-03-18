//!
//! Compression
//!
//! ## E-step
//!
//! Estimate node freq.
//!
//! ## M-step
//!
//! Solve min-flow and determine node copy number.
//!
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode};
use crate::hmmv2::freq::Reads;

///
///
///
pub fn compression<N: DbgNode, E: DbgEdge>(dbg: &Dbg<N, E>, reads: &Reads) -> Dbg<N, E> {
    // e-step
    // calculate node_freqs by using current dbg.

    // m-step
    // convert it to the
    unimplemented!();
}
