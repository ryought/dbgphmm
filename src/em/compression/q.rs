//!
//! Exact Q-function score
//!
use crate::common::CopyNum;
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::hmmv2::common::{PHMMEdge, PHMMNode};
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::{EdgeFreqs, NodeFreqs};

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
