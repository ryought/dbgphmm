//!
//! dbg as a seqgraph and phmm
//!
use super::dbg::{Dbg, DbgEdge, DbgEdgeBase, DbgNode, DbgNodeBase};
use crate::common::{CopyNum, Seq};
use crate::distribution::normal;
use crate::graph::seq_graph::{SeqEdge, SeqGraph, SeqNode};
use crate::hmmv2::common::PModel;
use crate::hmmv2::freq::PHMMOutput;
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::result::PHMMResultLike;
use crate::hmmv2::sample::State;
use crate::prob::Prob;
use rayon::prelude::*;

impl<N: DbgNode> SeqNode for N {
    fn copy_num(&self) -> CopyNum {
        self.copy_num()
    }
    fn base(&self) -> u8 {
        self.emission()
    }
}

impl<E: DbgEdge> SeqEdge for E {
    fn copy_num(&self) -> Option<CopyNum> {
        self.copy_num()
    }
}

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    /// Convert dbg into phmm
    pub fn to_phmm(&self, param: PHMMParams) -> PModel {
        self.graph.to_phmm(param)
    }
    ///
    /// Convert dbg into phmm and calculate full probability
    /// (Likelihood)
    ///
    pub fn to_full_prob<T>(&self, param: PHMMParams, seqs: T) -> Prob
    where
        T: IntoParallelIterator,
        T::Item: Seq,
    {
        let phmm = self.to_phmm(param);
        phmm.to_full_prob_parallel(seqs)
    }
    ///
    ///
    ///
    pub fn to_prior_prob(&self, genome_size_expected: CopyNum, genome_size_sigma: CopyNum) -> Prob {
        normal(
            self.genome_size() as f64,
            genome_size_expected as f64,
            genome_size_sigma as f64,
        )
    }
    pub fn to_prior_prob_by_lambda(&self, lambda: f64, genome_size_expected: CopyNum) -> Prob {
        let size_diff = genome_size_expected as f64 - self.genome_size() as f64;
        Prob::from_log_prob(-lambda * size_diff.powi(2))
    }
    ///
    /// calculate the prior score
    ///
    /// log P(G) = - lambda |G - G0|^2
    ///
    /// Parameters
    /// * lambda: penalty weight
    /// * G0: expected genome size
    ///
    pub fn to_prior_score_by_lambda(&self, lambda: f64, genome_size_expected: CopyNum) -> f64 {
        self.to_prior_prob_by_lambda(lambda, genome_size_expected)
            .to_log_value()
    }
}

impl<N: DbgNodeBase, E: DbgEdgeBase> Dbg<N, E> {
    ///
    /// show mapping of a sequence of emissions in cli.
    /// list the nodes which emits emissons[i] with high probability.
    ///
    pub fn show_mapping_summary<R: PHMMResultLike>(
        &self,
        emissions: &[u8],
        output: &PHMMOutput<R>,
    ) {
        let k = self.k();
        for (i, state_probs) in output.iter_emit_probs().skip(1).enumerate() {
            println!("{}{}", spaces(k + 1), emissions.to_str());
            println!("{}*{}", spaces(k + 1 + i), i);
            for (state, p) in state_probs.to_states().into_iter().take(10) {
                let ep = p.to_value();
                let s_state = match state {
                    State::Match(v) => format!("M {}", self.kmer(v)),
                    State::Ins(v) => format!("I {}", self.kmer(v)),
                    State::Del(v) => format!("D {}", self.kmer(v)),
                    State::MatchBegin => format!("MB"),
                    State::InsBegin => format!("IB"),
                    State::End => format!("E"),
                };
                let s_p = format!("{}", ep);
                println!("{}{} {}", spaces(i), s_state, s_p);
            }
        }
    }
}

///
/// get strings with repeated n-times space (' ').
///
fn spaces(n: usize) -> String {
    // old rust
    // std::iter::repeat(" ").take(n).collect::<String>()
    // new rust 1.16
    " ".repeat(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks::*;
    #[test]
    fn dbg_as_phmm_simple() {
        let dbg = mock_simple();
        let c = dbg.graph.total_emittable_copy_num();
        println!("c={}", c);
        let phmm = dbg.to_phmm(PHMMParams::default());
        println!("{}", phmm);
    }
}
