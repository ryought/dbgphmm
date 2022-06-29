//!
//! Exact Q-function score
//!
use crate::common::CopyNum;
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::hmmv2::common::{PHMMEdge, PHMMNode};
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::{EdgeFreqs, NodeFreqs};
use crate::min_flow::utils::clamped_log;

#[derive(Clone, Debug, Copy, Default)]
pub struct QScore {
    init: f64,
    trans: f64,
    prior: f64,
}

impl QScore {
    ///
    /// total score
    ///
    pub fn total(&self) -> f64 {
        self.init + self.trans + self.prior
    }
}

impl std::fmt::Display for QScore {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}(init={} trans={} prior={})",
            self.total(),
            self.init,
            self.trans,
            self.prior
        )
    }
}

///
/// calculate exact Q-function score with clamped log
///
pub fn q_score_clamped<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
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
                    * (clamped_log(copy_num) - clamped_log(copy_num_intersection));

                // (2) init
                qs.init += init_freq * (clamped_log(copy_num) - clamped_log(copy_num_total));
            }
        }
    }

    // (3) prior
    // -lambda (genome_size - genome_size_expected)^2
    let size_diff = genome_size as f64 - dbg.genome_size() as f64;
    qs.prior = -penalty_weight * size_diff.powi(2);

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
    let size_diff = genome_size as f64 - dbg.genome_size() as f64;
    qs.prior = -penalty_weight * size_diff.powi(2);

    qs
}
