//!
//! Functions to create multiple reads at once.
//!

///
/// Specifying read amount
///
#[derive(Clone, Debug)]
pub enum ReadAmount {
    ///
    /// By the number of reads.
    ///
    Count(usize),
    ///
    /// By the total bases of reads.
    ///
    TotalBases(usize),
}

impl PHMM {
    ///
    /// Generate a one-shot sequence of Emission and Hidden states
    /// by running a profile HMM.
    /// Random number generator will be created from the seed.
    ///
    pub fn sample_many(&self, length: usize, seed: u64, n_samples: usize) -> Historys {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        Historys(
            (0..n_samples)
                .map(|_| self.sample_rng(&mut rng, TerminateCond::Endable(length)))
                .collect(),
        )
    }
    ///
    /// Generate a sequence of emissions
    /// by running a profile HMM using rng(random number generator).
    ///
    /// (sampling reads from the model)
    ///
    pub fn sample_read(&self, length: usize, seed: u64) -> Sequence {
        let history = self.sample(length, seed);
        history.to_sequence()
    }
    ///
    /// Generate reads from sampling from profile
    ///
    pub fn sample_by_profile(&self, profile: &SampleProfile) -> Historys {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(profile.seed);
        let mut sampler_from = || match &profile.start_points {
            StartPoints::Custom(nodes) => {
                self.sample_rng_from_nodes(&mut rng, profile.length, nodes)
            }
            StartPoints::Random => self.sample_rng(&mut rng, profile.length),
            StartPoints::AllStartPoints => {
                panic!("StartPoints::AllStartPoints is not resolved until sample_by_profile")
            }
        };
        let mut sampler = || match profile.length {
            TerminateCond::EmitCount(count) => loop {
                let history = sampler_from();
                if history.total_bases() == count {
                    return history;
                }
            },
            _ => {
                let history = sampler_from();
                return history;
            }
        };
        let historys = match profile.read_amount {
            ReadAmount::Count(n_reads) => (0..n_reads).map(|_| sampler()).collect(),
            ReadAmount::TotalBases(required_total_bases) => {
                let mut historys = Vec::new();
                let mut total_bases = 0;
                while total_bases < required_total_bases {
                    let history = sampler();
                    total_bases += history.total_bases();
                    historys.push(history);
                }
                historys
            }
        };
        Historys(historys)
    }
}
