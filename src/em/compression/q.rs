//!
//! Exact Q-function score
//!
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase, NodeCopyNums};
use crate::dbg::float::FloatDbg;
use crate::graph::float_seq_graph::FloatSeqGraph;
use crate::hmmv2::common::{PHMMEdge, PHMMNode};
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::{EdgeFreqs, NodeFreqs};
use crate::min_flow::utils::clamped_log_with;
use crate::prelude::*;

// QScore definition was moved to hmmv2
pub use crate::hmmv2::q::QScore;

///
/// calculate exact Q-function score with clamped log
///
pub fn q_score_clamped<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
    zero_penalty: f64,
) -> QScore {
    let mut qs = QScore::default();

    let copy_num_total = dbg.genome_size();

    for intersection in dbg.iter_flow_intersections(edge_freqs) {
        let copy_num_intersection = intersection.total_in_copy_num();
        let freq_intersection = intersection.total_freq().unwrap();

        for node_info in intersection.iter_out_nodes() {
            let node = node_info.index;
            let copy_num = node_info.copy_num;
            let is_emittable = dbg.node(node).is_emittable();

            if is_emittable {
                let init_freq = init_freqs[node];
                let edge_freq_parents: f64 = dbg.parents(node).map(|(e, _, _)| edge_freqs[e]).sum();

                // calculate scores for the node
                // (1) trans
                qs.trans += edge_freq_parents
                    * (clamped_log_with(copy_num, zero_penalty)
                        - clamped_log_with(copy_num_intersection, zero_penalty));

                // (2) init
                qs.init += init_freq
                    * (clamped_log_with(copy_num, zero_penalty)
                        - clamped_log_with(copy_num_total, zero_penalty));
            }
        }
    }

    // (3) prior
    // -lambda (genome_size - genome_size_expected)^2
    qs.prior = dbg.to_prior_score(penalty_weight, genome_size);

    qs
}

///
/// Calculate (exact) Q function score.
///
pub fn q_score<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
) -> QScore {
    let mut qs = QScore::default();

    // the q_score does not depends on the phmmparams, so use the default one.
    let phmm = dbg.to_phmm(PHMMParams::default());

    for (node, node_weight) in phmm.nodes() {
        if phmm.is_emittable(node) {
            // (1) init score
            // A(Begin, v) log p_init(v)
            qs.init += init_freqs[node] * node_weight.init_prob().to_log_value();

            for (edge, child, edge_weight) in phmm.childs(node) {
                if phmm.is_emittable(child) {
                    // (2) trans score
                    // A(v,w) log p_trans(v, w)
                    qs.trans += edge_freqs[edge] * edge_weight.trans_prob().to_log_value();
                }
            }
        }
    }

    // (3) prior score
    // -lambda (genome_size - genome_size_expected)^2
    qs.prior = dbg.to_prior_score(penalty_weight, genome_size);

    qs
}

///
/// Calculate (exact) Q function score.
///
pub fn q_score_float<K: KmerLike>(
    dbg: &FloatDbg<K>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
) -> QScore {
    let mut init = 0.0;
    let mut trans = 0.0;

    for (node, node_weight) in dbg.nodes() {
        if node_weight.is_emittable() {
            // (1) init score
            // A(Begin, v) log p_init(v)
            init += init_freqs[node] * dbg.graph.init_prob(node).to_log_value();

            for (edge, child, edge_weight) in dbg.childs(node) {
                if dbg.is_emittable(child) {
                    // (2) trans score
                    // A(v,w) log p_trans(v, w)
                    trans += edge_freqs[edge] * dbg.graph.trans_prob(edge).to_log_value();
                }
            }
        }
    }

    // prior is always zero
    QScore::new(init, trans, 0.0)
}
