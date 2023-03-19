//!
//! Posterior probability inference of copy numbers on MultiDbg
//!
use super::{CopyNums, MultiDbg};
use crate::common::{CopyNum, PositionedReads, PositionedSequence, Seq};
use crate::distribution::normal;
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;
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
        unimplemented!();
    }
    ///
    ///
    ///
    pub fn max(&self) -> &CopyNums {
        unimplemented!();
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
}

///
/// Calculated score of copy numbers
///
/// Constructed by `MultiDbg::to_score`.
///
#[derive(Clone, Copy, Debug)]
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
    pub fn to_likelihood(&self, params: PHMMParams) -> Prob {
        unimplemented!();
    }
    ///
    /// Calculate the score `P(R|G)P(G)` (by prior `P(G)` from genome size and likelihood `P(R|G)` from reads) of this MultiDbg.
    ///
    pub fn to_score(
        &self,
        params: PHMMParams,
        genome_size_expected: CopyNum,
        genome_size_sigma: CopyNum,
    ) -> Score {
        Score {
            likelihood: self.to_likelihood(params),
            prior: self.to_prior(genome_size_expected, genome_size_sigma),
            genome_size: self.genome_size(),
            time: 0, // TODO
        }
    }
}

impl MultiDbg {
    ///
    /// # Arguments
    ///
    /// For likelihood
    /// * reads
    /// * PHMMParams
    ///
    /// For prior
    /// * genome_size_expected
    /// * genome_size_sigma
    ///
    /// For neighbors
    /// * max_cycle_size
    /// * max_flip
    ///
    /// # Procedure
    ///
    /// * evaluate scores of neighboring copy nums
    /// * move to the highest copy num
    /// * terminate if all neighbors have lower score
    ///
    pub fn sample_posterior(&self) -> Posterior {
        unimplemented!();
    }
}
