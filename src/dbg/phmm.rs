//!
//! dbg as a seqgraph and phmm
//!
use super::dbg::{Dbg, DbgEdge, DbgEdgeBase, DbgNode, DbgNodeBase};
use crate::common::{CopyNum, PositionedReads, PositionedSequence, Seq};
use crate::dbg::dbg::NodeCopyNums;
use crate::distribution::normal;
use crate::graph::seq_graph::{SeqEdge, SeqGraph, SeqNode};
use crate::hmmv2::common::PModel;
use crate::hmmv2::freq::PHMMOutput;
use crate::hmmv2::hint::Hint;
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::result::{PHMMResult, PHMMResultLike};
use crate::hmmv2::sample::State;
use crate::prob::Prob;
use crate::utils::{spaces, timer};
use itertools::Itertools;
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
    ///
    pub fn to_full_prob_sparse<T>(&self, param: PHMMParams, seqs: T) -> Prob
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
        let (p_rg, time) = timer(|| self.to_full_prob_sparse(param, seqs));
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

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
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
                state_probs.to_summary_string(|v| format!(
                    "{}x{}",
                    self.kmer(v).to_string(),
                    self.copy_num(v)
                ))
            );
            let p = output.forward.table_merged(i + 1).e().to_log_value();
            let p_prev = output.forward.table_merged(i).e().to_log_value();
            let dp = p - p_prev;
            println!("{}{}({})", spaces(i + 2), p, dp);
        }
    }
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
    pub fn compare_mappings_v2(
        &self,
        param: PHMMParams,
        emissions: &[u8],
        copy_nums_a: &NodeCopyNums,
        copy_nums_b: &NodeCopyNums,
    ) {
        let mut dbg_a = self.clone();
        dbg_a.set_node_copy_nums(copy_nums_a);
        let phmm_a = dbg_a.to_phmm(param);

        let mut dbg_b = self.clone();
        dbg_b.set_node_copy_nums(copy_nums_b);
        let phmm_b = dbg_b.to_phmm(param);

        let summary = |output: &PHMMOutput<PHMMResult>, copy_nums: &NodeCopyNums, x: usize| {
            let node_info = |v| format!("{}x{}", self.kmer(v).to_string(), copy_nums[v]);
            output.forward.tables[x]
                .to_states()
                .into_iter()
                .take(5)
                .map(|(state, prob)| {
                    if let Some(node) = state.to_node_index() {
                        format!("{}:{}{:.5}", node_info(node), state, prob.to_log_value())
                    } else {
                        format!("{}{:.5}", state, prob.to_log_value())
                    }
                })
                .join(",")
        };
        let dp = |output: &PHMMOutput<PHMMResult>, i: usize| {
            let p = output.forward.table_merged(i + 1).e().to_log_value();
            let p_prev = output.forward.table_merged(i).e().to_log_value();
            p - p_prev
        };

        let output_a = phmm_a.run(emissions);
        let output_b = phmm_b.run(emissions);
        let paf = output_a.to_full_prob_forward();
        let pab = output_a.to_full_prob_backward();
        let pbf = output_b.to_full_prob_forward();
        let pbb = output_b.to_full_prob_backward();
        println!("paf={} pab={} pbf={} pbb={}", paf, pab, pbf, pbb);

        let state_probs_a: Vec<_> = output_a.iter_emit_probs().skip(1).collect();
        let state_probs_b: Vec<_> = output_b.iter_emit_probs().skip(1).collect();

        let k = self.k();
        for i in 0..emissions.len() {
            // position
            let dp_a = dp(&output_a, i);
            let dp_b = dp(&output_b, i);
            // print header
            println!("{}{}", spaces(k + 7), emissions.to_str());
            println!(
                "{}i={:<5} {} {}",
                spaces(i),
                i,
                state_probs_a[i].to_summary_string_n(3, |v| format!(
                    "{}x{}x{}",
                    self.kmer(v).to_string(),
                    copy_nums_a[v],
                    copy_nums_b[v],
                )),
                dp_a,
            );
            println!(
                "{}i={:<5} {} {}",
                spaces(i),
                i,
                state_probs_b[i].to_summary_string_n(3, |v| format!(
                    "{}x{}x{}",
                    self.kmer(v).to_string(),
                    copy_nums_a[v],
                    copy_nums_b[v],
                )),
                dp_b,
            );
            println!("{} {:.5}", spaces(i), dp_b - dp_a);
        }
    }
    pub fn compare_mappings(
        &self,
        param: PHMMParams,
        reads: &PositionedReads,
        copy_nums_a: &NodeCopyNums,
        copy_nums_b: &NodeCopyNums,
    ) {
        let mut dbg_a = self.clone();
        dbg_a.set_node_copy_nums(copy_nums_a);
        let phmm_a = dbg_a.to_phmm(param);

        let mut dbg_b = self.clone();
        dbg_b.set_node_copy_nums(copy_nums_b);
        let phmm_b = dbg_b.to_phmm(param);

        let summary = |output: &PHMMOutput<PHMMResult>, copy_nums: &NodeCopyNums, x: usize| {
            let node_info = |v| format!("{}x{}", self.kmer(v).to_string(), copy_nums[v]);
            output.forward.tables[x]
                .to_states()
                .into_iter()
                .take(5)
                .map(|(state, prob)| {
                    if let Some(node) = state.to_node_index() {
                        format!("{}:{}{:.5}", node_info(node), state, prob.to_log_value())
                    } else {
                        format!("{}{:.5}", state, prob.to_log_value())
                    }
                })
                .join(",")
        };

        for (i, read) in reads.into_iter().enumerate() {
            let emissions = read.seq();
            let output_a = phmm_a.run(emissions);
            let output_b = phmm_b.run(emissions);
            let pa = output_a.to_full_prob_forward().to_log_value();
            let pb = output_b.to_full_prob_forward().to_log_value();
            println!(
                "r[{}]\tsummary\t{}\t{}\t{}\t{}",
                i,
                read.origin_node().index(),
                pa,
                pb,
                pa - pb,
            );

            for x in 0..emissions.len() {
                let pa = output_a.forward.table_merged(x + 1).e().to_log_value();
                let pa_prev = output_a.forward.table_merged(x).e().to_log_value();
                let dpa = pa - pa_prev;
                let sa = summary(&output_a, &copy_nums_a, x);

                let pb = output_b.forward.table_merged(x + 1).e().to_log_value();
                let pb_prev = output_b.forward.table_merged(x).e().to_log_value();
                let dpb = pb - pb_prev;
                let sb = summary(&output_b, &copy_nums_b, x);

                println!(
                    "r[{}]\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    i,
                    x,
                    emissions[x] as char,
                    pa,
                    pb,
                    dpa,
                    dpb,
                    dpa - dpb,
                    sa,
                    sb
                );
            }
            // self.show_mapping_summary(read.as_ref(), &output);
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
