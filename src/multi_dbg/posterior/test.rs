//!
//! Test of posterior sampling
//!
use super::super::MultiDbg;
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
    use crate::hmmv2::params::PHMMParams;
    use crate::utils::timer;

    #[test]
    #[ignore]
    fn simple_genome() {
        let (genome, genome_size) = genome::simple(1000, 0);
        let param = PHMMParams::uniform(0.005);
        let dataset = generate_dataset(
            genome.clone(),
            genome_size,
            0,
            20,
            500,
            ReadType::FragmentWithRevComp,
            param,
        );
        dataset.show_reads_with_genome();

        let k = 20;
        let (mdbg, t) = timer(|| MultiDbg::create_draft_from_dataset(k, &dataset));
        mdbg.to_gfa_file("simple.gfa");
        mdbg.to_dbg_file("simple.dbg");

        let paths_true = mdbg.compact_paths_from_styled_seqs(&genome).unwrap();
        mdbg.to_paths_file("simple.paths", &paths_true);

        let param_infer = PHMMParams::uniform(0.01);
        let post = mdbg.sample_posterior_with_dataset(&dataset, param_infer, 200, 10, 1, 10);
        post.to_file("simple.post");
        mdbg.to_gfa_post_file("simple.post.gfa", &post);

        // compare with true copynums
    }

    #[test]
    #[ignore]
    fn repeat_1k() {
        let (genome, genome_size) =
            genome::tandem_repeat_polyploid_with_unique_homo_ends(100, 10, 0, 300, 2, 0.01, 0);
        let param = PHMMParams::uniform(0.01);
        let dataset = generate_dataset(
            genome.clone(),
            genome_size,
            0,
            20,
            500,
            ReadType::FragmentWithRevComp,
            param,
        );
        dataset.show_reads_with_genome();

        let k = 20;
        let (mdbg, t) = timer(|| MultiDbg::create_draft_from_dataset(k, &dataset));
        println!("created mdbg in {}ms", t);
        mdbg.to_gfa_file("reads.gfa");
        mdbg.to_dbg_file("reads.dbg");

        let param_infer = PHMMParams::uniform(0.01);
        mdbg.sample_posterior_with_dataset(&dataset, param_infer, 200, 10, 1, 10);
    }
}
