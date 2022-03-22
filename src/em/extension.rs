//!
//! Extension
//!
//! ## E-step
//!
//!
//! ## M-step
//!
//!
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode};
use crate::hmmv2::freq::Reads;
use crate::hmmv2::params::PHMMParams;

///
/// Extension algorithm
///
/// ## Details
///
/// * e-step
///     Calculate edge_freqs on dbg.
///
/// * m-step
///     Maximize the score for each intersections
///
/// ## TODOs
///
/// * avoid dbg copy?
///
pub fn extension<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
) -> Dbg<N, E> {
    unimplemented!();
}
