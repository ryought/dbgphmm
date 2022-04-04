//!
//! **Extension**
//!
//! convert k-dbg into k+1-dbg.
//!
//! ## E-step
//!
//!
//! ## M-step
//!
//!
use crate::common::{CopyNum, Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, EdgeCopyNums};
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::trans_table::EdgeFreqs;
use crate::kmer::kmer::KmerLike;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
pub mod flow_intersection;
use flow_intersection::{FlowIntersection, FlowIntersectionEdge, FlowIntersectionNode};
pub mod intersection_graph;

///
/// Extension algorithm
///
/// ## Details
///
/// * e-step
///     Calculate edge_freqs on dbg.
///
/// * m-step
///     Maximize the score for each intersections.
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
    // (1) infer edge freqs
    // (2) infer the best copy nums
    unimplemented!();
}

///
/// E-step of extension
///
/// ## Details
///
/// * convert dbg into phmm.
/// * calculate edge frequencies by forward/backward algorithm to emit the reads.
///
fn extension_e_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
) -> EdgeFreqs {
    let phmm = dbg.to_phmm(params.clone());
    phmm.to_edge_freqs(reads)
}

///
/// M-step of extension
///
/// ## Params
///
/// * `dbg`
///     de bruijn graph
/// * `edge_freqs`
///     edge usage frequencies.
///
/// ## Details
///
///
fn extension_m_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
) -> EdgeCopyNums {
    unimplemented!();
}
