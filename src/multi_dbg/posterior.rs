//!
//! Posterior probability inference of copy numbers on MultiDbg
//!
use super::{
    neighbors::{NeighborConfig, UpdateInfo, UpdateMethod},
    CopyNums, MultiDbg, Path,
};
use crate::common::{CopyNum, PositionedReads, PositionedSequence, ReadCollection, Seq};
use crate::distribution::normal;
use crate::e2e::Dataset;
use crate::hist::DiscreteDistribution;
use crate::hmmv2::hint::Mappings;
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;
use crate::utils::{progress_common_style, timer};
use fnv::FnvHashMap as HashMap;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use petgraph::graph::{EdgeIndex, NodeIndex};
use rayon::prelude::*;
use rustflow::min_flow::residue::{
    update_cycle_from_str, update_cycle_to_string, ResidueDirection, UpdateCycle,
};
use serde::{Deserialize, Serialize};

pub mod output;
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
    pub samples: Vec<PosteriorSample>,
    ///
    /// Total probability of sampled copy numbers
    ///
    pub p: Prob,
}

///
///
///
#[derive(Clone, Debug, PartialEq)]
pub struct PosteriorSample {
    pub copy_nums: CopyNums,
    pub score: Score,
    pub infos: Vec<UpdateInfo>,
}

impl PosteriorSample {
    ///
    ///
    ///
    fn to_infos_string_internal(infos: &[UpdateInfo]) -> String {
        format!("[{}]", infos.iter().format(","))
    }
    ///
    ///
    ///
    pub fn to_infos_string(&self) -> String {
        Self::to_infos_string_internal(&self.infos)
    }
    ///
    /// Convert `[e5+e2-e4+,e5+e6+]` into Vec<UpdateCycle>
    ///
    pub fn from_infos_str(s: &str) -> Option<Vec<UpdateInfo>> {
        let s = s.strip_prefix('[')?;
        let s = s.strip_suffix(']')?;
        let mut infos = Vec::new();
        for e in s.split(',') {
            if !e.is_empty() {
                let info = e.parse().ok()?;
                infos.push(info);
            }
        }
        Some(infos)
    }
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
    pub fn add(&mut self, sample: PosteriorSample) {
        if !self.contains(&sample.copy_nums) {
            self.p += sample.score.p();
            self.samples.push(sample);
        }
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
    pub fn find(&self, copy_nums: &CopyNums) -> Option<&PosteriorSample> {
        self.samples
            .iter()
            .find(|sample| &sample.copy_nums == copy_nums)
    }
    ///
    /// Get the best posterior sample with highest score.
    ///
    pub fn max_sample(&self) -> &PosteriorSample {
        self.samples
            .iter()
            .max_by_key(|sample| sample.score.p())
            .unwrap()
    }
    ///
    /// Get the best copy numbers with highest score.
    ///
    pub fn max_copy_nums(&self) -> &CopyNums {
        let sample = self.max_sample();
        &sample.copy_nums
    }
    ///
    ///
    ///
    pub fn samples(&self) -> &[PosteriorSample] {
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
            .map(|sample| (sample.copy_nums[edge], sample.score.p() / self.p()))
            .collect();
        DiscreteDistribution::from_occurs(&copy_nums_with_prob)
    }
}

///
/// Calculated score of copy numbers
///
/// Constructed by `MultiDbg::to_score`.
///
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Score {
    ///
    /// Likelihood `P(R|X)`
    ///
    pub likelihood: Prob,
    ///
    /// Prior `P(G)`
    ///
    pub prior: Prob,
    ///
    /// Genome size `|G|`
    ///
    pub genome_size: CopyNum,
    ///
    /// Number of Euler circuits of DBG (stored in log)
    ///
    /// `P(X) = P(G) * n_euler_circuits`
    ///
    pub n_euler_circuits: f64,
    ///
    /// Computation time of likelihood (alignment)
    ///
    pub time_likelihood: u128,
    ///
    /// Computation time of n_euler_circuit (mainly matrix determinant)
    ///
    pub time_euler: u128,
}

impl Score {
    ///
    /// Calculate total probability `P(R|X)P(X)`
    ///
    pub fn p(&self) -> Prob {
        self.likelihood * self.prior * Prob::from_log_prob(self.n_euler_circuits)
    }
}

impl std::fmt::Display for Score {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", serde_json::to_string(self).unwrap())
    }
}

impl std::str::FromStr for Score {
    type Err = serde_json::Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        serde_json::from_str(s)
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
    pub fn to_likelihood<S: Seq>(
        &self,
        param: PHMMParams,
        reads: &ReadCollection<S>,
        mappings: Option<&Mappings>,
    ) -> Prob {
        let phmm = self.to_phmm(param);
        phmm.to_full_prob_reads(reads, mappings, true)
    }
    ///
    /// Calculate the score `P(R|G)P(G)` (by prior `P(G)` from genome size and likelihood `P(R|G)` from reads) of this MultiDbg.
    ///
    pub fn to_score<S: Seq>(
        &self,
        param: PHMMParams,
        reads: &ReadCollection<S>,
        mappings: Option<&Mappings>,
        genome_size_expected: CopyNum,
        genome_size_sigma: CopyNum,
    ) -> Score {
        let (likelihood, time_likelihood) = timer(|| self.to_likelihood(param, reads, mappings));
        let (n_euler_circuits, time_euler) = timer(|| self.n_euler_circuits());
        Score {
            likelihood,
            prior: self.to_prior(genome_size_expected, genome_size_sigma),
            genome_size: self.genome_size(),
            n_euler_circuits,
            time_likelihood,
            time_euler,
        }
    }
}

///
/// Posterior sampling function
///
/// * Mapping: `generate_hints`
/// * Sampling: `sample_posterior`
/// * Upgrading: `purge_and_extend_with_posterior`
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
        mappings: &Mappings,
        genome_size_expected: CopyNum,
        genome_size_sigma: CopyNum,
        neighbor_config: NeighborConfig,
        max_iter: usize,
        rescue_only: bool,
    ) -> Posterior {
        let mut post = Posterior::new();
        let mut copy_nums = self.get_copy_nums();
        let mut infos = Vec::new();
        let mut dbg = self.clone();
        let mut n_iter = 0;
        let freqs = mappings.to_node_freqs(dbg.n_edges_full());
        let coverage = reads.total_bases() as f64 / genome_size_expected as f64;
        println!("sample_posterior rescue_only={}", rescue_only);

        // calculate initial score
        let score = dbg.to_score(
            param,
            reads,
            Some(mappings),
            genome_size_expected,
            genome_size_sigma,
        );
        post.add(PosteriorSample {
            copy_nums: copy_nums.clone(),
            score,
            infos: Vec::new(),
        });

        'outer: while n_iter < max_iter {
            println!("n_iter={}", n_iter);
            // calculate scores of new neighboring copynums of current copynum
            //
            dbg.set_copy_nums(&copy_nums);

            // A. Generate neighbors
            //
            // [0] rescue 0x -> 1x changes
            let rescue_neighbors = dbg.to_rescue_neighbors(&freqs, coverage, 2, 10, true);

            let neighbor_copy_nums_set = if rescue_only {
                vec![rescue_neighbors]
            } else {
                // [1] partial search
                let partial_neighbors = dbg.to_neighbor_copy_nums_and_infos(NeighborConfig {
                    max_cycle_size: 5,
                    max_flip: 2,
                    use_long_cycles: true,
                    ignore_cycles_passing_terminal: true,
                    use_reducers: false,
                });
                // [2] full search
                let full_neighbors = dbg.to_neighbor_copy_nums_and_infos(neighbor_config);
                vec![rescue_neighbors, partial_neighbors, full_neighbors]
            };

            // B. Try each neighbor and move to the better neighbor if found.
            //
            for (i, neighbor_copy_nums) in neighbor_copy_nums_set.into_iter().enumerate() {
                eprintln!(
                    "iter #{}-set{} n_neighbors={}",
                    n_iter,
                    i,
                    neighbor_copy_nums.len()
                );
                match dbg.sample_posterior_once(
                    neighbor_copy_nums,
                    &mut post,
                    &infos,
                    param,
                    reads,
                    mappings,
                    genome_size_expected,
                    genome_size_sigma,
                    rescue_only, // enable multi-move mode if rescue only
                ) {
                    Some(sample) => {
                        eprintln!("iter #{}-set{} early terminate", n_iter, i);
                        eprintln!("accepted={}", sample.to_infos_string());
                        copy_nums = sample.copy_nums;
                        infos = sample.infos;
                        n_iter += 1;
                        continue 'outer;
                    }
                    None => {}
                }
            }

            // C. if no better neighbor could not be found (i.e. current copynums is local optimum)
            // terminate sampling.
            //
            eprintln!("iter #{} not found", n_iter);
            break 'outer;
        }

        post
    }
    ///
    ///
    ///
    pub fn sample_posterior_for_inspect<S: Seq>(
        &self,
        neighbors: Vec<(CopyNums, UpdateInfo)>,
        param: PHMMParams,
        reads: &ReadCollection<S>,
        genome_size_expected: CopyNum,
        genome_size_sigma: CopyNum,
    ) -> Posterior {
        // evaluate (calculate) score of the given copy numbers and return a PosteriorSample
        let evaluate = |copy_nums: CopyNums, infos: Vec<UpdateInfo>| {
            // evaluate score
            let mut dbg = self.clone();
            dbg.set_copy_nums(&copy_nums);
            let score = dbg.to_score(param, reads, None, genome_size_expected, genome_size_sigma);
            PosteriorSample {
                copy_nums,
                score,
                infos,
            }
        };

        let mut posterior = Posterior::new();

        // init
        posterior.add(evaluate(self.get_copy_nums(), vec![]));

        // neighbors
        let samples: Vec<_> = neighbors
            .into_par_iter()
            .progress_with_style(progress_common_style())
            .filter_map(|(copy_nums, info)| {
                eprintln!("info={:?}", info);
                if posterior.contains(&copy_nums) {
                    None
                } else {
                    Some(evaluate(copy_nums, vec![info]))
                }
            })
            .collect();

        for sample in samples.into_iter() {
            posterior.add(sample);
        }

        posterior
    }
    ///
    ///
    ///
    pub fn sample_posterior_once<S: Seq>(
        &self,
        neighbors: Vec<(CopyNums, UpdateInfo)>,
        posterior: &mut Posterior,
        infos_init: &[UpdateInfo],
        param: PHMMParams,
        reads: &ReadCollection<S>,
        mappings: &Mappings,
        genome_size_expected: CopyNum,
        genome_size_sigma: CopyNum,
        multi_move: bool,
    ) -> Option<PosteriorSample> {
        // evaluate (calculate) score of the given copy numbers and return a PosteriorSample
        let evaluate = |copy_nums: CopyNums, info: UpdateInfo| {
            // evaluate score
            let mut dbg = self.clone();
            dbg.set_copy_nums(&copy_nums);
            let score = dbg.to_score(
                param,
                reads,
                Some(mappings),
                genome_size_expected,
                genome_size_sigma,
            );
            let mut infos = infos_init.to_owned();
            infos.push(info);
            PosteriorSample {
                copy_nums,
                score,
                infos,
            }
        };

        let t_start = std::time::Instant::now();
        let samples: Vec<_> = neighbors
            .clone()
            .into_par_iter()
            .progress_with_style(progress_common_style())
            .filter_map(|(copy_nums, info)| {
                if posterior.contains(&copy_nums) {
                    None
                } else {
                    Some(evaluate(copy_nums, info))
                }
            })
            .collect();
        let time = t_start.elapsed().as_millis();
        println!(
            "sampled n_samples={} k={} n_reads={} total_bases={} in t={}ms",
            samples.len(),
            self.k(),
            reads.len(),
            reads.total_bases(),
            time
        );
        for sample in samples.iter() {
            posterior.add(sample.clone());
        }

        //
        // accept the best neighbor
        //
        if multi_move {
            // accept independent moves
            let mut current_copy_nums = self.get_copy_nums();
            let original_copy_nums = self.get_copy_nums();
            let mut accepted_update_cycles = vec![];
            let current_score = posterior
                .find(&current_copy_nums)
                .expect("current copy number was not sampled")
                .score;
            eprintln!("multi move from p={}", current_score.p());

            for (i, (copy_nums, info)) in neighbors
                .iter()
                .sorted_by_key(|(copy_nums, _)| posterior.find(&copy_nums).unwrap().score.p())
                .rev()
                .enumerate()
            {
                let score = posterior.find(&copy_nums).unwrap().score;

                eprintln!(
                    "{} {} {} {}",
                    i,
                    update_cycle_to_string(&info.cycle()),
                    score.p(),
                    info.cycle()
                        .iter()
                        .map(|(edge, _)| format!(
                            "e{}={}x",
                            edge.index(),
                            original_copy_nums[*edge]
                        ))
                        .join(",")
                );

                // if the change improves the score and independent (does not conflicts with other
                // accepted changes)
                let improves_score = score.p() > current_score.p();
                let cycle = info.cycle();
                let is_independent =
                    self.is_independent_update(&current_copy_nums, &accepted_update_cycles, &cycle);
                if improves_score && is_independent {
                    eprintln!("accept!! {} {}", improves_score, is_independent);
                    self.apply_update_info_to_copy_nums(&mut current_copy_nums, &cycle);
                    accepted_update_cycles.push(cycle.clone());
                } else {
                    eprintln!("reject {} {}", improves_score, is_independent);
                }

                // neighbors are sorted by score
                if !improves_score {
                    break;
                }
            }

            // run once and check if the score actually improves.
            let info = UpdateInfo::new(accepted_update_cycles, UpdateMethod::MultiMove);
            let sample = evaluate(current_copy_nums, info);
            println!("multi move score {}", sample.score.p());
            posterior.add(sample);
        }

        // move to highest copy num and continue
        // if better copynums was found, move to it.
        let sample = posterior.max_sample();
        if sample.copy_nums != self.get_copy_nums() {
            Some(sample.clone())
        } else {
            None
        }
    }
    ///
    /// Append hint information for reads in parallel
    ///
    /// # Procedure
    ///
    /// 1. Convert MultiDbg into PHMM with uniform transition probability
    /// 2. Run forward/backward on the PHMM and obtain node freqs of each read
    ///
    pub fn generate_mappings<S: Seq>(
        &self,
        mut param: PHMMParams,
        reads: &ReadCollection<S>,
        mappings: Option<&Mappings>,
    ) -> Mappings {
        param.n_warmup = self.k();
        // let phmm = self.to_uniform_phmm(param);
        let phmm = self.to_non_zero_phmm(param);
        let (map, time) = timer(|| phmm.generate_mappings(reads, mappings, true, Some(50.0)));
        eprintln!(
            "generated mappings for k={} n_reads={} total_bases={} in t={}ms",
            self.k(),
            reads.len(),
            reads.total_bases(),
            time
        );

        // self.inspect_freqs(&map, 20.0);

        map
    }
    pub fn inspect_freqs(&self, mappings: &Mappings, coverage: f64) {
        let freqs = self.mappings_to_freqs(&mappings);
        let copy_num = self.min_squared_error_copy_nums_from_freqs(&freqs, coverage, None);
        for edge_compact in self.graph_compact().edge_indices() {
            println!(
                "e{} {}x {}x",
                edge_compact.index(),
                self.copy_num_of_edge_in_compact(edge_compact),
                copy_num[edge_compact]
            );
            for (i, &edge_full) in self.edges_in_full(edge_compact).iter().enumerate() {
                println!("\t{} {}", i, freqs[NodeIndex::new(edge_full.index())]);
            }
        }
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
    pub fn purge_and_extend_with_posterior(
        &self,
        posterior: &Posterior,
        k_max: usize,
        p0: Prob,
        paths: Option<Vec<Path>>,
        mappings: &Mappings,
    ) -> (Self, Option<Vec<Path>>, Mappings) {
        // (0)
        // Set copy number to the most probable one
        let mut dbg = self.clone();
        dbg.set_copy_nums(posterior.max_copy_nums());

        // (1)
        // Find edges to be purged according to posterior distribution
        // List edges whose current copynum is 0 and posterior probability P(X=0) is high
        //
        let mut edges_purge = Vec::new();
        for edge in dbg.graph_compact().edge_indices() {
            if posterior.p_edge_x(edge, 0) > p0 && dbg.copy_num_of_edge_in_compact(edge) == 0 {
                edges_purge.push(edge);
            }
        }

        // (2)
        // Do purge and extend
        let (dbg, paths, mappings) =
            dbg.purge_and_extend(&edges_purge, k_max, true, paths, mappings);
        (dbg, paths, mappings)
    }
}

///
/// Get posterior distribution of DBG of k_max
///
/// 1. posterior for current k-DBG
/// 2. purge 0x edges
/// 3. extend to k+1
///
pub fn infer_posterior_by_extension<
    S: Seq,
    F: Fn(&MultiDbg, &Posterior, &Option<Vec<Path>>, &Mappings),
    G: Fn(&MultiDbg, &Mappings),
    H: Fn(&MultiDbg, &Mappings),
>(
    k_max: usize,
    dbg_init: MultiDbg,
    // evaluate
    param_infer: PHMMParams,
    param_error: PHMMParams,
    reads: &ReadCollection<S>,
    // prior
    genome_size_expected: CopyNum,
    genome_size_sigma: CopyNum,
    // neighbor
    neighbor_config: NeighborConfig,
    max_iter: usize,
    // extend
    p0: Prob,
    // callback
    on_iter: F,
    on_map: G,
    on_extend: H,
    // true path if available
    paths: Option<Vec<Path>>,
    // inference parameters
    k_max_rescue_only: usize,
    k_max_rerun_mapping: usize,
    // precalculated mapping
    mappings: Option<Mappings>,
) -> (MultiDbg, Posterior, Option<Vec<Path>>, Mappings) {
    let mut dbg = dbg_init;
    let mut mappings = match mappings {
        Some(mappings) => {
            eprintln!("using precalculated mappings");
            mappings
        }
        None => {
            eprintln!("generating mappings");
            dbg.generate_mappings(param_error, reads, None)
        }
    };
    on_map(&dbg, &mappings);
    let mut paths = paths;
    let mut posterior;
    let coverage = reads.total_bases() as f64 / genome_size_expected as f64;

    loop {
        eprintln!("k={}", dbg.k());

        // (2) posterior
        let t_start_posterior = std::time::Instant::now();
        posterior = dbg.sample_posterior(
            param_infer,
            &reads,
            &mappings,
            genome_size_expected,
            genome_size_sigma,
            neighbor_config,
            max_iter,
            dbg.k() <= k_max_rescue_only,
        );
        dbg.set_copy_nums(posterior.max_copy_nums());
        let t_posterior = t_start_posterior.elapsed();
        eprintln!(
            "posterior sampling k={} t={}ms",
            dbg.k(),
            t_posterior.as_millis()
        );

        // (3) run callback
        on_iter(&dbg, &posterior, &paths, &mappings);

        if dbg.k() >= k_max {
            // dbg does not need extension
            break;
        }

        // (4) extend
        let t_start_extend = std::time::Instant::now();
        (dbg, paths, mappings) =
            dbg.purge_and_extend_with_posterior(&posterior, k_max, p0, paths, &mappings);
        let t_extend = t_start_extend.elapsed();
        eprintln!("extend t={}ms", t_extend.as_millis());
        on_extend(&dbg, &mappings);

        // (1) update hints before extending
        let t_start_hint = std::time::Instant::now();
        if dbg.k() <= k_max_rerun_mapping {
            println!("k={} rerun mapping", dbg.k());
            mappings = dbg.generate_mappings(param_error, reads, None); // currently previous mapping
                                                                        // is not used
        } else {
            println!("k={} not rerun mapping", dbg.k());
            mappings = dbg.generate_mappings(param_error, reads, Some(&mappings));
        }
        let t_hint = t_start_hint.elapsed();
        eprintln!("hint t={}ms", t_hint.as_millis());
        on_map(&dbg, &mappings);

        // (1b) approximate copy numbers from the mapping and frequency
        let t_start_approx = std::time::Instant::now();
        let freqs = dbg.mappings_to_freqs(&mappings);
        let copy_num = dbg.min_squared_error_copy_nums_from_freqs(&freqs, coverage, Some(2));
        dbg.set_copy_nums(&copy_num);
        let t_approx = t_start_approx.elapsed();
        eprintln!("approx t={}ms coverage={}", t_approx.as_millis(), coverage);
    }

    // final run using p_error
    let t_start_posterior = std::time::Instant::now();
    let mappings = dbg.generate_mappings(param_error, reads, None);
    posterior = dbg.sample_posterior(
        param_error,
        &reads,
        &mappings,
        genome_size_expected,
        genome_size_sigma,
        neighbor_config,
        max_iter,
        false,
    );
    dbg.set_copy_nums(posterior.max_copy_nums());
    let t_posterior = t_start_posterior.elapsed();
    eprintln!("posterior_final t={}ms", t_posterior.as_millis());

    (dbg, posterior, paths, mappings)
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ei;
    use crate::prob::p;

    #[test]
    fn score() {
        let a = Score {
            likelihood: Prob::from_prob(0.3),
            prior: Prob::from_prob(0.5),
            genome_size: 111,
            time_likelihood: 102,
            n_euler_circuits: 10.0,
            time_euler: 100,
        };
        let t = a.to_string();
        println!("{}", t);
        let b: Score = t.parse().unwrap();
        println!("{}", b);
        assert_eq!(a, b);
    }
}
