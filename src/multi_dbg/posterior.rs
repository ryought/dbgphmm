//!
//! Posterior probability inference of copy numbers on MultiDbg
//!
use super::{CopyNums, MultiDbg};
use crate::common::{CopyNum, PositionedReads, PositionedSequence, ReadCollection, Seq};
use crate::distribution::normal;
use crate::hist::DiscreteDistribution;
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;
use crate::utils::timer;
use petgraph::graph::EdgeIndex;

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
    pub fn p_edge(&self, edge: EdgeIndex, x: CopyNum) -> Prob {
        unimplemented!();
    }
    ///
    /// Posterior distribution of copy number of the edge `P(X[edge] | R)`
    ///
    pub fn p_edge_dist(&self, edge: EdgeIndex) -> DiscreteDistribution {
        unimplemented!();
    }
}

///
/// dump and load functions
///
/// ```text
/// P   -19281.0228
/// C   -192919.0    likelihood=0.0    [1,2,1,1,1,2,1,0]
/// C   -193799.0    likelihood=0.0    [1,2,1,1,0,2,0,0]
/// ```
///
impl Posterior {
    ///
    ///
    ///
    pub fn to_writer<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        writeln!(writer, "P\t{}", self.p.to_log_value())?;
        for (copy_nums, score) in self.samples.iter() {
            writeln!(
                writer,
                "C\t{}\t{}\t{}",
                score.p().to_log_value(),
                score,
                copy_nums
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
                    let score: Score = iter.next().unwrap().parse().unwrap();
                    let copy_nums: CopyNums = iter.next().unwrap().parse().unwrap();

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
    ) -> Posterior {
        let mut post = Posterior::new();
        let mut copy_nums = self.get_copy_nums();
        let mut dbg = self.clone();
        let mut n_iter = 0;

        while n_iter < max_iter {
            // calculate scores of new neighboring copynums of current copynum
            //
            dbg.set_copy_nums(&copy_nums);
            for (copy_nums, _info) in dbg.to_neighbor_copy_nums_and_infos(max_cycle_size, max_flip)
            {
                if !post.contains(&copy_nums) {
                    // evaluate score
                    dbg.set_copy_nums(&copy_nums);
                    let score = dbg.to_score(param, reads, genome_size_expected, genome_size_sigma);
                    post.add(copy_nums, score);
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
                likelihood: p(0.001),
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
