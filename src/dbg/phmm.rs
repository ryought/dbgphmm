//!
//! dbg as a seqgraph and phmm
//!
use super::dbg::{Dbg, DbgEdge, DbgEdgeBase, DbgNode, DbgNodeBase};
use crate::common::{CopyNum, Seq};
use crate::distribution::normal;
use crate::graph::seq_graph::{SeqEdge, SeqGraph, SeqNode};
use crate::hmmv2::common::PModel;
use crate::hmmv2::freq::PHMMOutput;
use crate::hmmv2::hint::Hint;
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::result::PHMMResultLike;
use crate::hmmv2::sample::State;
use crate::prob::Prob;
use crate::utils::{spaces, timer};
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

//
// EvalResult
//
#[derive(Clone, Copy)]
pub struct EvalResult {
    /// Likelihood `P(R|G)`
    p_rg: Prob,
    /// Prior `P(G)`
    p_g: Prob,
    /// Genome size `|G|`
    genome_size: CopyNum,
    /// Computation time of likelihood
    time: u128,
}
impl EvalResult {
    ///
    /// Constructor of EvalResult
    ///
    /// (used in Dbg::evaluate)
    ///
    fn new(p_rg: Prob, p_g: Prob, genome_size: CopyNum, time: u128) -> Self {
        EvalResult {
            p_rg,
            p_g,
            genome_size,
            time,
        }
    }
    ///
    /// Unnormalized posterior P(R,G) = P(R|G)P(G)
    ///
    pub fn posterior(&self) -> Prob {
        self.p_rg * self.p_g
    }
    /// compuation time
    pub fn time(&self) -> u128 {
        self.time
    }
    /// Likelihood P(R|G)
    pub fn p_rg(&self) -> Prob {
        self.p_rg
    }
    /// Prior P(G)
    pub fn p_g(&self) -> Prob {
        self.p_g
    }
}
impl std::fmt::Display for EvalResult {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "P(R|G)={} P(G={})={} T={}",
            self.p_rg, self.genome_size, self.p_g, self.time,
        )
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
    /// Convert dbg into phmm and calculate full probability
    /// (Likelihood)
    ///
    pub fn to_full_prob_with_hint<S>(&self, param: PHMMParams, seqs_and_hints: &[(S, Hint)]) -> Prob
    where
        S: Seq,
    {
        let phmm = self.to_phmm(param);
        phmm.to_full_prob_par_with_hint(seqs_and_hints)
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
    ///
    ///
    ///
    pub fn evaluate<T>(
        &self,
        param: PHMMParams,
        seqs: T,
        genome_size_expected: CopyNum,
        genome_size_sigma: CopyNum,
    ) -> EvalResult
    where
        T: IntoParallelIterator,
        T::Item: Seq,
    {
        let (p_rg, time) = timer(|| self.to_full_prob(param, seqs));
        let p_g = self.to_prior_prob(genome_size_expected, genome_size_sigma);
        EvalResult::new(p_rg, p_g, self.genome_size(), time)
    }
    ///
    ///
    ///
    pub fn evaluate_with_hint<S: Seq>(
        &self,
        param: PHMMParams,
        seqs_and_hints: &[(S, Hint)],
        genome_size_expected: CopyNum,
        genome_size_sigma: CopyNum,
    ) -> EvalResult {
        let (p_rg, time) = timer(|| self.to_full_prob_with_hint(param, seqs_and_hints));
        let p_g = self.to_prior_prob(genome_size_expected, genome_size_sigma);
        EvalResult::new(p_rg, p_g, self.genome_size(), time)
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
        let header = || {
            println!("{}{}", spaces(k + 7), emissions.to_str());
        };
        for (i, state_probs) in output.iter_emit_probs().skip(1).enumerate() {
            // print header
            if i % 10 == 0 {
                header();
            }
            println!(
                "{}i={:<5} {}",
                spaces(i),
                i,
                state_probs.to_summary_string(|v| self.kmer(v).to_string())
            );
            let p = output.forward.table_merged(i + 1).e().to_log_value();
            let p_prev = output.forward.table_merged(i).e().to_log_value();
            let dp = p - p_prev;
            println!("{}{}({})", spaces(i + 2), p, dp);
        }
    }
}
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    pub fn show_mapping_summary_for_reads<T>(&self, param: PHMMParams, reads: T)
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        let phmm = self.to_phmm(param);
        for (i, read) in reads.into_iter().enumerate() {
            let output = phmm.run(read.as_ref());
            println!("r[{}] {}", i, read.as_ref().to_str());
            self.show_mapping_summary(read.as_ref(), &output);
        }
    }
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
    #[test]
    fn dbg_phmm_show_mapping_summary() {
        let dbg = mock_simple();
        let read = b"CTTGCTT";
        let phmm = dbg.to_phmm(PHMMParams::default());
        let output = phmm.run(read);
        dbg.show_mapping_summary(read, &output);
    }
}
