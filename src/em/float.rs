//!
//! Float dbg optimization by EM & GradDescent
//!
//! * takes FloatDbg as input
//! * E-step: run forward/backward algorithm and calculate the node/edge usage
//! * M-step: improve Q score by GradDescent and MinCostFlow
//! * iterate E/M-steps and returns the improved FloatDbg
//!
use crate::dbg::float::{CopyDensity, FloatDbg, FloatDbgEdge, FloatDbgNode};
use crate::graph::float_seq_graph::FloatSeqGraph;
use crate::hmmv2::{EdgeFreqs, NodeFreqs};
use crate::prelude::*;

///
/// E-step: calculate edge_freqs (freq between v->w) and init_freqs (freq between Begin->w)
///
pub fn e_step<K: KmerLike>(
    dbg: &FloatDbg<K>,
    reads: &Reads,
    params: &PHMMParams,
) -> (EdgeFreqs, NodeFreqs, Prob) {
    let phmm = dbg.graph.to_phmm(params.clone());
    phmm.to_edge_and_init_freqs_parallel(reads)
}

///
/// M-step
///
pub fn m_step<K: KmerLike>(
    dbg: &FloatDbg<K>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyDensity,
) {
    // (1) convert to edge-centric dbg with each edge has a cost

    // (2) search for negative cycle

    // search again for
}

///
/// FloatDbg
/// -> construct residue graph of QScore difference when +1/-1
/// -> change density if negative meaningful cycle was found
/// -> new updated FloatDbg
///
pub fn m_step_once<K: KmerLike>(
    dbg: &FloatDbg<K>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyDensity,
) -> FloatDbg<K> {
    unimplemented!();
}

fn to_residue_graph() {}

struct ScoreDiffs {
    /// difference amount
    diff: f64,
    /// difference of score when flow is increased by diff
    diff_inc: f64,
    /// difference of score when flow is decreased by diff
    diff_dec: f64,
}