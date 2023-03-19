//!
//! Posterior probability inference of copy numbers on MultiDbg
//!
use super::{CopyNums, MultiDbg};
use crate::common::{CopyNum, PositionedReads, PositionedSequence, Seq};
use crate::distribution::normal;
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;

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

///
///
///
#[derive(Clone, Copy, Debug)]
pub struct Score {
    /// Likelihood `P(R|G)`
    likelihood: Prob,
    /// Prior `P(G)`
    prior: Prob,
    /// Genome size `|G|`
    genome_size: CopyNum,
    /// Computation time of likelihood
    time: u128,
}

///
/// Scoreing related functions
///
/// * to_likelihood
/// * to_prior
/// * to_score
///
impl MultiDbg {
    /// Calculate the prior probability `P(G)`
    ///
    /// genome size follows Normal distribution with mean = genome_size_expected and var =
    /// genome_size_sigma.
    ///
    pub fn to_prior_prob(&self, genome_size_expected: CopyNum, genome_size_sigma: CopyNum) -> Prob {
        normal(
            self.genome_size() as f64,
            genome_size_expected as f64,
            genome_size_sigma as f64,
        )
    }
    /// Calculate the likelihood `P(R|G)`
    ///
    /// convert MultiDbg into Profile HMM (PHMMModel) and calculate the full probability of reads.
    ///
    pub fn to_likelihood(&self, params: PHMMParams) -> Prob {
        unimplemented!();
    }
    /// Calculate the score (prior and likelihood) of this MultiDbg.
    ///
    pub fn to_score(&self) -> Score {
        unimplemented!();
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
    ///
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
