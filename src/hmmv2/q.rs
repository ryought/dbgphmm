//!
//! Q score calculation functions
//! after forward/backward calculation
//!
//! - q_score
//! - q_score_diff
//!
use crate::hmmv2::common::{PHMMEdge, PHMMModel, PHMMNode};
use crate::hmmv2::{EdgeFreqs, NodeFreqs};

#[derive(Clone, Debug, Copy, Default)]
pub struct QScore {
    /// init score
    pub init: f64,
    /// trans score
    pub trans: f64,
    /// prior score
    pub prior: f64,
}

impl QScore {
    ///
    /// Constructor
    ///
    pub fn new(init: f64, trans: f64, prior: f64) -> Self {
        QScore { init, trans, prior }
    }
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
/// Calculate (exact) Q function score.
///
pub fn q_score_exact<N: PHMMNode, E: PHMMEdge>(
    phmm: &PHMMModel<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
) -> QScore {
    let mut init = 0.0;
    let mut trans = 0.0;

    for (node, node_weight) in phmm.nodes() {
        if node_weight.is_emittable() {
            // (1) init score
            // A(Begin, v) log p_init(v)
            let init_prob = node_weight.init_prob().to_log_value();
            assert!(init_prob.is_finite());
            init += init_freqs[node] * init_prob;

            // (2) trans score
            // A(v,w) log p_trans(v, w)
            for (edge, child, edge_weight) in phmm.childs(node) {
                if phmm.is_emittable(child) {
                    let trans_prob = edge_weight.trans_prob().to_log_value();
                    assert!(trans_prob.is_finite());
                    trans += edge_freqs[edge] * trans_prob;
                }
            }
        }
    }

    // prior is always zero
    QScore::new(init, trans, 0.0)
}
