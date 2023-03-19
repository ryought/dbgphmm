//!
//! Test of posterior sampling
//!
use crate::e2e::Dataset;

///
/// 1. Create draft MultiDbg
/// 2. Infer posterior distribution
/// 3. Check the posterior probability of true copynums
/// 4. Check the posterior probability of copy num of edge of false-kmer
///
fn test_posterior(dataset: &Dataset) {}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::e2e::{generate_dataset, ReadType};
    use crate::genome;

    #[test]
    fn simple_genome() {
        genome::simple(1000, 0);
    }
}
