//!
//! Posterior probability inference of copy numbers on MultiDbg
//!
use super::{CopyNums, MultiDbg, Path};
use crate::common::{CopyNum, PositionedReads, PositionedSequence, ReadCollection, Seq};
use crate::distribution::normal;
use crate::e2e::Dataset;
use crate::hist::DiscreteDistribution;
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;
use crate::utils::timer;
use fnv::FnvHashMap as HashMap;
use indicatif::{ParallelProgressIterator, ProgressStyle};
use itertools::Itertools;
use petgraph::graph::EdgeIndex;
use rayon::prelude::*;

pub mod test;

///
/// Collection of sampled CopyNums and its scores.
///
/// * samples
/// * p
///
#[derive(Clone, Debug)]
pub struct Posterior {
    ///
    /// Collection of sampled copy nums and its score
    ///
    samples: Vec<(CopyNums, Score)>,
    ///
    /// Total probability of sampled copy numbers
    ///
    p: Prob,
}

impl Posterior {
    ///
    /// Create empty posterior container
    ///
    pub fn new() -> Self {
        Posterior {
            samples: Vec::new(),
            p: Prob::zero(),
        }
    }
    ///
    /// Add a sampled copy numbers and its score
    ///
    pub fn add(&mut self, copy_nums: CopyNums, score: Score) {
        self.samples.push((copy_nums, score));
        self.p += score.p();
    }
    ///
    /// Check if the copy numbers is stored in the posterior or not
    ///
    pub fn contains(&self, copy_nums: &CopyNums) -> bool {
        self.find(copy_nums).is_some()
    }
    ///
    ///
    ///
    pub fn find(&self, copy_nums: &CopyNums) -> Option<&(CopyNums, Score)> {
        self.samples
            .iter()
            .find(|(copy_nums_in, _)| copy_nums_in == copy_nums)
    }
    ///
    /// Get the best copy numbers with highest score.
    ///
    pub fn max_copy_nums(&self) -> &CopyNums {
        let (copy_nums, _score) = self
            .samples
            .iter()
            .max_by_key(|(_copy_nums, score)| score.p())
            .unwrap();
        copy_nums
    }
    ///
    ///
    ///
    pub fn samples(&self) -> &[(CopyNums, Score)] {
        &self.samples
    }
    ///
    /// Sum of total probability: normalization factor of posterior probability
    ///
    pub fn p(&self) -> Prob {
        self.p
    }
    ///
    /// Posterior probability of copy number of the edge
    ///
    /// `P(X[edge] = x | R)`
    ///
    pub fn p_edge_x(&self, edge: EdgeIndex, x: CopyNum) -> Prob {
        self.p_edge(edge).p_x(x)
    }
    ///
    /// Posterior distribution of copy number of the edge `P(X[edge] | R)`
    ///
    pub fn p_edge(&self, edge: EdgeIndex) -> DiscreteDistribution {
        let copy_nums_with_prob: Vec<_> = self
            .samples
            .iter()
            .map(|(copy_nums, score)| (copy_nums[edge], score.p() / self.p()))
            .collect();
        DiscreteDistribution::from_occurs(&copy_nums_with_prob)
    }
}

///
/// dump and load functions
///
/// ```text
/// Z   -19281.0228
/// C   -192919.0    [1,2,1,1,1,2,1,0]   likelihood=0.00
/// C   -191882.0    [1,2,0,0,1,2,2,1]   likelihood=0.01
/// ```
///
impl Posterior {
    ///
    ///
    ///
    pub fn to_writer<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        writeln!(writer, "Z\t{}", self.p.to_log_value())?;
        for (copy_nums, score) in self
            .samples
            .iter()
            .sorted_by_key(|(_, score)| score.p())
            .rev()
        {
            writeln!(
                writer,
                "C\t{}\t{}\t{}",
                score.p().to_log_value(),
                copy_nums,
                score,
            )?
        }
        Ok(())
    }
    ///
    /// create string with `to_gfa_writer`
    ///
    pub fn to_string(&self) -> String {
        let mut writer = Vec::with_capacity(128);
        self.to_writer(&mut writer).unwrap();
        String::from_utf8(writer).unwrap()
    }
    ///
    /// create file with `to_gfa_writer`
    ///
    pub fn to_file<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        let mut file = std::fs::File::create(path).unwrap();
        self.to_writer(&mut file)
    }
    ///
    ///
    ///
    pub fn from_reader<R: std::io::BufRead>(reader: R) -> Self {
        let mut samples = Vec::new();
        let mut p = Prob::zero();

        for line in reader.lines() {
            let text = line.unwrap();
            let first_char = text.chars().nth(0).unwrap();
            match first_char {
                'C' => {
                    let mut iter = text.split_whitespace();
                    iter.next().unwrap(); // 'C'
                    iter.next().unwrap(); // value
                    let copy_nums: CopyNums = iter.next().unwrap().parse().unwrap();
                    let score: Score = iter.next().unwrap().parse().unwrap();

                    samples.push((copy_nums, score));
                    p += score.p();
                }
                _ => {} // ignore
            }
        }

        Posterior { samples, p }
    }
    ///
    ///
    pub fn from_str(s: &str) -> Self {
        Self::from_reader(s.as_bytes())
    }
    ///
    ///
    pub fn from_file<P: AsRef<std::path::Path>>(path: P) -> Self {
        let file = std::fs::File::open(path).unwrap();
        let reader = std::io::BufReader::new(file);
        Self::from_reader(reader)
    }
}

///
/// benchmark functions for when true genome is available
///
impl MultiDbg {
    ///
    /// Everytime
    /// * posterior probability (normalized)
    /// * likelihood (log)
    /// * prior (log)
    /// * genome size
    ///
    /// Only if genome is known
    /// * diff of copynums from true
    ///
    pub fn to_inspect_writer<W: std::io::Write>(
        &self,
        mut writer: W,
        posterior: &Posterior,
        copy_nums_true: &CopyNums,
    ) -> std::io::Result<()> {
        // for each copy nums
        for (i, (copy_nums, score)) in posterior
            .samples
            .iter()
            .sorted_by_key(|(_, score)| score.p())
            .rev()
            .enumerate()
        {
            writeln!(
                writer,
                "{}\tC\t{}\t{:.10}\t{}\t{}\t{}\t{}\t{}",
                self.k(),
                i,
                (score.p() / posterior.p()).to_value(),
                score.likelihood.to_log_value(),
                score.prior.to_log_value(),
                score.genome_size,
                copy_nums.diff(&copy_nums_true),
                copy_nums,
            )?
        }

        // for each edges
        for edge in self.graph_compact().edge_indices() {
            let p_edge = posterior.p_edge(edge);
            let copy_num_true = copy_nums_true[edge];
            writeln!(
                writer,
                "{}\tE\te{}\t{}\t{:.5}\t{:.5}\t{:.5}\t{}",
                self.k(),
                edge.index(),
                copy_num_true,
                p_edge.mean(),
                p_edge.p_x(copy_num_true).to_value(),
                p_edge.p_x(0).to_value(),
                p_edge.to_short_string(),
            )?
        }

        Ok(())
    }
    ///
    ///
    pub fn to_inspect_string(&self, posterior: &Posterior, copy_nums_true: &CopyNums) -> String {
        let mut writer = Vec::with_capacity(128);
        self.to_inspect_writer(&mut writer, posterior, copy_nums_true)
            .unwrap();
        String::from_utf8(writer).unwrap()
    }
    ///
    ///
    pub fn to_inspect_file<P: AsRef<std::path::Path>>(
        &self,
        path: P,
        posterior: &Posterior,
        copy_nums_true: &CopyNums,
    ) -> std::io::Result<()> {
        let mut file = std::fs::File::create(path).unwrap();
        self.to_inspect_writer(&mut file, posterior, copy_nums_true)
    }
}

///
/// Calculated score of copy numbers
///
/// Constructed by `MultiDbg::to_score`.
///
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Score {
    ///
    /// Likelihood `P(R|G)`
    ///
    likelihood: Prob,
    ///
    /// Prior `P(G)`
    ///
    prior: Prob,
    ///
    /// Genome size `|G|`
    ///
    genome_size: CopyNum,
    ///
    /// Computation time of likelihood
    ///
    time: u128,
}

impl Score {
    ///
    /// Calculate total probability `P(R|G)P(G)`
    ///
    fn p(&self) -> Prob {
        self.likelihood * self.prior
    }
}

impl std::fmt::Display for Score {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "likelihood={},prior={},genome_size={},time={}",
            self.likelihood, self.prior, self.genome_size, self.time
        )
    }
}

impl std::str::FromStr for Score {
    type Err = std::num::ParseFloatError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut likelihood = None;
        let mut prior = None;
        let mut genome_size = None;
        let mut time = None;

        for e in s.split(',') {
            let mut it = e.split('=');
            let key = it.next().unwrap();
            let value = it.next().unwrap();
            match key {
                "likelihood" => {
                    likelihood = Some(value.parse().unwrap());
                }
                "prior" => {
                    prior = Some(value.parse().unwrap());
                }
                "genome_size" => {
                    genome_size = Some(value.parse().unwrap());
                }
                "time" => {
                    time = Some(value.parse().unwrap());
                }
                _ => {}
            }
        }

        Ok(Score {
            likelihood: likelihood.unwrap(),
            prior: prior.unwrap(),
            genome_size: genome_size.unwrap(),
            time: time.unwrap(),
        })
    }
}

///
/// Scoreing related functions
///
/// * to_likelihood
/// * to_prior
/// * to_score
///
impl MultiDbg {
    ///
    /// Calculate the prior probability `P(G)`
    ///
    /// genome size follows Normal distribution with mean = genome_size_expected and var =
    /// genome_size_sigma.
    ///
    pub fn to_prior(&self, genome_size_expected: CopyNum, genome_size_sigma: CopyNum) -> Prob {
        normal(
            self.genome_size() as f64,
            genome_size_expected as f64,
            genome_size_sigma as f64,
        )
    }
    ///
    /// Calculate the likelihood `P(R|G)`
    ///
    /// convert MultiDbg into Profile HMM (PHMMModel) and calculate the full probability of reads.
    ///
    pub fn to_likelihood<S: Seq>(&self, param: PHMMParams, reads: &ReadCollection<S>) -> Prob {
        let phmm = self.to_phmm(param);
        phmm.to_full_prob_reads(reads)
    }
    ///
    /// Calculate the score `P(R|G)P(G)` (by prior `P(G)` from genome size and likelihood `P(R|G)` from reads) of this MultiDbg.
    ///
    pub fn to_score<S: Seq>(
        &self,
        param: PHMMParams,
        reads: &ReadCollection<S>,
        genome_size_expected: CopyNum,
        genome_size_sigma: CopyNum,
    ) -> Score {
        let (likelihood, time) = timer(|| self.to_likelihood(param, reads));
        Score {
            likelihood,
            prior: self.to_prior(genome_size_expected, genome_size_sigma),
            genome_size: self.genome_size(),
            time,
        }
    }
}

///
/// Posterior sampling function
///
impl MultiDbg {
    ///
    /// # Arguments
    ///
    /// ## For likelihood
    /// * reads
    /// * PHMMParams
    ///
    /// ## For prior
    /// * genome_size_expected
    /// * genome_size_sigma
    ///
    /// ## For neighbors
    /// * max_cycle_size
    /// * max_flip
    ///
    /// ## For greedy search
    /// * max_iter
    ///
    ///
    /// # Procedure
    ///
    /// * start from current copynums
    /// * evaluate scores of neighboring copynums
    /// * move to the highest copynums
    /// * terminate if all neighbors have lower score
    ///
    pub fn sample_posterior<S: Seq>(
        &self,
        param: PHMMParams,
        reads: &ReadCollection<S>,
        genome_size_expected: CopyNum,
        genome_size_sigma: CopyNum,
        max_cycle_size: usize,
        max_flip: usize,
        max_iter: usize,
        is_parallel: bool,
    ) -> Posterior {
        let mut post = Posterior::new();
        let mut copy_nums = self.get_copy_nums();
        let mut dbg = self.clone();
        let mut n_iter = 0;

        while n_iter < max_iter {
            // calculate scores of new neighboring copynums of current copynum
            //
            dbg.set_copy_nums(&copy_nums);
            let neighbor_copy_nums = dbg.to_neighbor_copy_nums_and_infos(max_cycle_size, max_flip);
            eprintln!("iter#{} n_neighbors={}", n_iter, neighbor_copy_nums.len());

            // evaluate all neighbors
            if is_parallel {
                let style = ProgressStyle::with_template(
                    "[{elapsed_precise}/{eta_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
                )
                .unwrap()
                .progress_chars("##-");

                let neighbors_with_score: Vec<_> = neighbor_copy_nums
                    .into_par_iter()
                    .progress_with_style(style)
                    .filter_map(|(copy_nums, _info)| {
                        if post.contains(&copy_nums) {
                            None
                        } else {
                            // evaluate score
                            let mut dbg = self.clone();
                            dbg.set_copy_nums(&copy_nums);
                            let score =
                                dbg.to_score(param, reads, genome_size_expected, genome_size_sigma);
                            Some((copy_nums, score))
                        }
                    })
                    .collect();
                for (copy_nums, score) in neighbors_with_score {
                    post.add(copy_nums, score);
                }
            } else {
                for (i, (copy_nums, _info)) in neighbor_copy_nums.into_iter().enumerate() {
                    eprintln!("iter#{} neighbor#{}", n_iter, i);
                    if !post.contains(&copy_nums) {
                        // evaluate score
                        dbg.set_copy_nums(&copy_nums);
                        let score =
                            dbg.to_score(param, reads, genome_size_expected, genome_size_sigma);
                        post.add(copy_nums, score);
                    }
                }
            }

            // move to highest copy num and continue
            // if self is highest, terminate
            let next_copy_nums = post.max_copy_nums().clone();
            if next_copy_nums == copy_nums {
                break;
            } else {
                copy_nums = next_copy_nums;
                n_iter += 1;
            }
        }

        post
    }
    ///
    /// Append hint information for reads in parallel
    ///
    /// # Procedure
    ///
    /// 1. Convert MultiDbg into PHMM with uniform transition probability
    /// 2. Run forward/backward on the PHMM and obtain node freqs of each read
    ///
    pub fn generate_hints<S: Seq>(
        &self,
        param: PHMMParams,
        reads: ReadCollection<S>,
        parallel: bool,
    ) -> ReadCollection<S> {
        let phmm = self.to_uniform_phmm(param);
        phmm.append_hints(reads, parallel)
    }
    /// Extend to k+1 by sampled posterior distribution
    ///
    /// 1. Purge edges by posterior
    /// 2. Extend Dbg into k+1
    /// 3. Convert Hints and Path for k+1 if necessary
    ///
    /// # Arguments
    ///
    /// * `posterior`: sampled posterior distribution
    /// * `p0`: remove edge `e` if `P(X(e)=0|R) > p0`
    /// * `paths` (optional)
    /// * `reads` (optional)
    ///
    pub fn extend_with_posterior<S: Seq>(
        &self,
        posterior: &Posterior,
        k_max: usize,
        p0: Prob,
        paths: Option<Vec<Path>>,
        reads: ReadCollection<S>,
    ) -> (Self, Option<Vec<Path>>, ReadCollection<S>) {
        // (1)
        // Find edges to be purged according to posterior distribution
        // List edges whose current copynum is 0 and posterior probability P(X=0) is high
        //
        let mut edges_purge = Vec::new();
        for edge in self.graph_compact().edge_indices() {
            if posterior.p_edge_x(edge, 0) > p0 && self.copy_num_of_edge_in_compact(edge) == 0 {
                edges_purge.push(edge);
            }
        }

        // (2)
        // Do purge and extend
        if reads.has_hint() {
            let hints = reads.hints.unwrap();
            let (dbg_kp1, paths, hints) =
                self.purge_and_extend(&edges_purge, k_max, paths, Some(hints));
            let reads = ReadCollection::from_with_hint(reads.reads, hints.unwrap());
            (dbg_kp1, paths, reads)
        } else {
            let (dbg_kp1, paths, _) = self.purge_and_extend(&edges_purge, k_max, paths, None);
            let reads = ReadCollection::from(reads.reads);
            (dbg_kp1, paths, reads)
        }
    }
}

///
/// Get posterior distribution of DBG of k_max
///
/// 1. posterior for current k-DBG
/// 2. purge 0x edges
/// 3. extend to k+1
///
pub fn infer_posterior_by_extension<S: Seq, F: Fn(MultiDbg, Posterior, Option<Vec<Path>>)>(
    k_max: usize,
    dbg_init: MultiDbg,
    // evaluate
    param: PHMMParams,
    reads: ReadCollection<S>,
    genome_size_expected: CopyNum,
    genome_size_sigma: CopyNum,
    // neighbor
    max_cycle_size: usize,
    max_flip: usize,
    max_iter: usize,
    // extend
    p0: Prob,
    // callback
    on_extend: F,
    // true path if available
    paths: Option<Vec<Path>>,
) -> (MultiDbg, Posterior, Option<Vec<Path>>) {
    let mut dbg = dbg_init;
    let mut reads = reads;
    let mut paths = paths;

    /*
    let mut posterior;
    loop {
        posterior = dbg.sample_posterior(
            param,
            &reads,
            genome_size_expected,
            genome_size_sigma,
            max_cycle_size,
            max_flip,
            max_iter,
            true,
        );

        while dbg.k() < k_final {
            // ignore if unneccessary
            let is_ambiguous = dbg.n_ambiguous_node() > 0;
            (dbg, paths, reads) = dbg.extend_with_posterior(&posterior, p0, paths, reads);
            if is_ambiguous {
                break;
            }
        }
    }
    (dbg, posterior, paths)
    */

    unimplemented!();
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prob::p;

    #[test]
    fn posterior_dump_load() {
        let mut post = Posterior::new();
        post.add(
            vec![1, 1, 2, 2, 0].into(),
            Score {
                likelihood: p(0.6),
                prior: p(0.3),
                time: 10,
                genome_size: 101,
            },
        );
        post.add(
            vec![1, 1, 1, 2, 1].into(),
            Score {
                likelihood: p(0.003),
                prior: p(0.2),
                time: 11,
                genome_size: 99,
            },
        );
        let s = post.to_string();
        println!("{}", s);

        let post_loaded = Posterior::from_str(&s);
        assert_eq!(post_loaded.samples, post.samples);
        assert_eq!(post_loaded.p, post.p);
    }

    #[test]
    fn score() {
        let a = Score {
            likelihood: Prob::from_prob(0.3),
            prior: Prob::from_prob(0.5),
            genome_size: 111,
            time: 102,
        };
        let t = a.to_string();
        println!("{}", t);
        let b: Score = t.parse().unwrap();
        println!("{}", b);
        assert_eq!(a, b);
    }
}
