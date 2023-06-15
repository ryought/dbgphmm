//!
//! Posterior probability inference of copy numbers on MultiDbg
//!
use super::{neighbors::NeighborConfig, CopyNums, MultiDbg, Path};
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
use rustflow::min_flow::residue::{ResidueDirection, UpdateInfo};

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
    pub fn to_infos_string(&self) -> String {
        format!(
            "[{}]",
            self.infos
                .iter()
                .map(|info| {
                    info.iter()
                        .map(|(edge, dir)| format!("e{}{}", edge.index(), dir))
                        .join("")
                })
                .join(",")
        )
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
/// dump and load functions
///
/// ```text
/// Z   -19281.0228
/// C   -192919.0    [1,2,1,1,1,2,1,0]   likelihood=0.00 e1+e2-e3+,e1+e2-
/// C   -191882.0    [1,2,0,0,1,2,2,1]   likelihood=0.01
/// ```
///
impl Posterior {
    ///
    ///
    ///
    pub fn to_writer<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        writeln!(writer, "# {}", env!("GIT_HASH"))?;
        writeln!(writer, "Z\t{}", self.p.to_log_value())?;
        for sample in self
            .samples
            .iter()
            .sorted_by_key(|sample| sample.score.p())
            .rev()
        {
            writeln!(
                writer,
                "C\t{}\t{}\t{}\t{}",
                sample.score.p().to_log_value(),
                sample.copy_nums,
                sample.score,
                sample.to_infos_string(),
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
                    let infos: Vec<UpdateInfo> = iter
                        .next()
                        .unwrap()
                        .trim_start_matches('[')
                        .trim_end_matches(']')
                        .split_terminator(',')
                        .map(|s| {
                            let mut info = Vec::new();
                            for x in s.split_inclusive(&['+', '-']) {
                                let index: usize = x
                                    .trim_start_matches('e')
                                    .trim_end_matches(&['+', '-'])
                                    .parse()
                                    .unwrap();
                                let dir: ResidueDirection =
                                    x.rmatches(&['+', '-']).next().unwrap().parse().unwrap();
                                info.push((EdgeIndex::new(index), dir));
                            }
                            info
                        })
                        .collect();
                    p += score.p();
                    samples.push(PosteriorSample {
                        copy_nums,
                        score,
                        infos,
                    });
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

fn format_option_copy_num<T: std::fmt::Display>(o: Option<T>) -> String {
    match o {
        Some(i) => format!("{}", i),
        None => format!("?"),
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
        copy_nums_true: Option<&CopyNums>,
    ) -> std::io::Result<()> {
        // for each copy nums
        // TODO clean up using key value format
        writeln!(writer, "# {}", env!("GIT_HASH"))?;
        writeln!(
            writer,
            "{}\tG\tn_edges_full\t{}",
            self.k(),
            self.n_edges_full()
        )?;
        writeln!(
            writer,
            "{}\tG\tn_edges_compact\t{}",
            self.k(),
            self.n_edges_compact()
        )?;
        writeln!(
            writer,
            "{}\tG\tn_nodes_full\t{}",
            self.k(),
            self.n_nodes_full()
        )?;
        writeln!(
            writer,
            "{}\tG\tn_nodes_compact\t{}",
            self.k(),
            self.n_nodes_compact()
        )?;
        writeln!(
            writer,
            "{}\tG\tn_emittable_edges\t{}",
            self.k(),
            self.n_emittable_edges()
        )?;
        writeln!(
            writer,
            "{}\tG\tdegree_stats\t{:?}",
            self.k(),
            self.degree_stats(),
        )?;
        for (i, sample) in posterior
            .samples
            .iter()
            .sorted_by_key(|sample| sample.score.p())
            .rev()
            .enumerate()
        {
            let score = &sample.score;
            let copy_nums = &sample.copy_nums;
            writeln!(
                writer,
                "{}\tC\t{}\t{:.10}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.k(),
                i,
                (score.p() / posterior.p()).to_value(),
                score.likelihood.to_log_value(),
                score.prior.to_log_value(),
                score.n_euler_circuits,
                score.genome_size,
                format_option_copy_num(
                    copy_nums_true.map(|copy_nums_true| copy_nums_true.diff(copy_nums))
                ),
                sample.to_infos_string(),
                copy_nums,
            )?
        }

        // for each edges
        for edge in self.graph_compact().edge_indices() {
            let p_edge = posterior.p_edge(edge);
            let copy_num_true = copy_nums_true.map(|copy_num_true| copy_num_true[edge]);
            writeln!(
                writer,
                "{}\tE\te{}\t{}\t{:.5}\t{:.5}\t{:.5}\t{}",
                self.k(),
                edge.index(),
                format_option_copy_num(copy_num_true),
                p_edge.mean(),
                format_option_copy_num(
                    copy_num_true.map(|copy_num_true| p_edge.p_x(copy_num_true).to_value())
                ),
                p_edge.p_x(0).to_value(),
                p_edge.to_short_string(),
            )?
        }

        Ok(())
    }
    ///
    ///
    pub fn to_inspect_string(
        &self,
        posterior: &Posterior,
        copy_nums_true: Option<&CopyNums>,
    ) -> String {
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
        copy_nums_true: Option<&CopyNums>,
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
        write!(
            f,
            "likelihood={},prior={},genome_size={},n_euler_circuits={},time_likelihood={},time_euler={}",
            self.likelihood, self.prior, self.genome_size, self.n_euler_circuits, self.time_likelihood, self.time_euler
        )
    }
}

impl std::str::FromStr for Score {
    type Err = std::num::ParseFloatError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut likelihood = None;
        let mut prior = None;
        let mut genome_size = None;
        let mut n_euler_circuits = None;
        let mut time_likelihood = None;
        let mut time_euler = None;

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
                "n_euler_circuits" => {
                    n_euler_circuits = Some(value.parse().unwrap());
                }
                "time_likelihood" => {
                    time_likelihood = Some(value.parse().unwrap());
                }
                "time_euler" => {
                    time_euler = Some(value.parse().unwrap());
                }
                _ => {}
            }
        }

        Ok(Score {
            likelihood: likelihood.unwrap(),
            prior: prior.unwrap(),
            genome_size: genome_size.unwrap(),
            n_euler_circuits: n_euler_circuits.unwrap(),
            time_likelihood: time_likelihood.unwrap(),
            time_euler: time_euler.unwrap(),
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
    pub fn to_likelihood<S: Seq>(
        &self,
        param: PHMMParams,
        reads: &ReadCollection<S>,
        mappings: Option<&Mappings>,
    ) -> Prob {
        let phmm = self.to_phmm(param);
        phmm.to_full_prob_reads(reads, mappings)
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
            let rescue_neighbors = dbg.to_rescue_neighbors(2, 10);
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
                    false,
                    // rescue_only, // enable multi-move mode if rescue only
                ) {
                    Some(sample) => {
                        eprintln!("iter #{}-set{} early terminate", n_iter, i);
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
        let t_start = std::time::Instant::now();
        let samples: Vec<_> = neighbors
            .into_par_iter()
            .progress_with_style(progress_common_style())
            .filter_map(|(copy_nums, info)| {
                if posterior.contains(&copy_nums) {
                    None
                } else {
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
                    Some(PosteriorSample {
                        copy_nums,
                        score,
                        infos,
                    })
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
            unimplemented!();
        } else {
            // move to highest copy num and continue
            // if better copynums was found, move to it.
            let sample = posterior.max_sample();
            if sample.copy_nums != self.get_copy_nums() {
                Some(sample.clone())
            } else {
                None
            }
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
        let phmm = self.to_uniform_phmm(param);
        let (map, time) = timer(|| phmm.generate_mappings(reads, mappings, true));
        println!(
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
            dbg.generate_mappings(param_infer, reads, None)
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
            mappings = dbg.generate_mappings(param_infer, reads, None); // currently previous mapping
                                                                        // is not used
        } else {
            println!("k={} not rerun mapping", dbg.k());
            mappings = dbg.generate_mappings(param_infer, reads, Some(&mappings));
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
    fn posterior_dump_load() {
        let mut post = Posterior::new();
        post.add(PosteriorSample {
            copy_nums: vec![1, 1, 2, 2, 0].into(),
            score: Score {
                likelihood: p(0.6),
                prior: p(0.3),
                time_likelihood: 10,
                genome_size: 101,
                n_euler_circuits: 10.0,
                time_euler: 12,
            },
            infos: vec![
                vec![
                    (ei(0), ResidueDirection::Up),
                    (ei(121), ResidueDirection::Down),
                ],
                vec![
                    (ei(3), ResidueDirection::Down),
                    (ei(1), ResidueDirection::Up),
                ],
            ],
        });
        post.add(PosteriorSample {
            copy_nums: vec![1, 1, 1, 2, 1].into(),
            score: Score {
                likelihood: p(0.003),
                prior: p(0.2),
                time_likelihood: 11,
                genome_size: 99,
                n_euler_circuits: 10.0,
                time_euler: 12,
            },
            infos: Vec::new(),
        });
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
