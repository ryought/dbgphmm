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
use crate::common::{Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode};
use crate::hmmv2::params::PHMMParams;
use crate::min_flow::Cost;
use crate::prob::Prob;
pub mod v1;
pub mod v2;
pub mod v3;
pub use v1::compression_step;

///
/// Log information store of each iteration in compression
///
#[derive(Clone)]
pub struct CompressionLog<N: DbgNode, E: DbgEdge> {
    /// Full probability
    pub full_prob: Prob,
    /// Min-flow error
    pub min_flow_score: Cost,
    /// resulting dbg
    pub dbg: Dbg<N, E>,
}

impl<N: DbgNode, E: DbgEdge> CompressionLog<N, E> {
    pub fn new(full_prob: Prob, min_flow_score: Cost, dbg: Dbg<N, E>) -> Self {
        CompressionLog {
            full_prob,
            min_flow_score,
            dbg,
        }
    }
}

impl<N: DbgNode, E: DbgEdge> std::fmt::Display for CompressionLog<N, E> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}",
            self.full_prob.to_log_value(),
            self.min_flow_score,
            self.dbg
        )
    }
}

///
/// Compression full algorithm by running `compression_step` iteratively.
///
/// * max_iter: max iteration loop count of EM.
///
pub fn compression<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    depth: Freq,
    max_iter: usize,
) -> (Dbg<N, E>, Vec<CompressionLog<N, E>>) {
    let depths = vec![depth; max_iter];
    compression_with_depths(dbg, reads, params, &depths)
}

///
/// Compression full algorithm by running `compression_step` iteratively
/// with specified depths.
///
pub fn compression_with_depths<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    depths: &[Freq],
) -> (Dbg<N, E>, Vec<CompressionLog<N, E>>) {
    let mut dbg = dbg.clone();
    let mut logs = Vec::new();

    // iterate EM steps
    for &depth in depths {
        let (dbg_new, is_updated, log) = v1::compression_step(&dbg, reads, params, depth);
        logs.push(log);

        // if the single EM step does not change the DBG model, stop iteration.
        if !is_updated {
            break;
        }
        dbg = dbg_new;
    }

    (dbg, logs)
}
