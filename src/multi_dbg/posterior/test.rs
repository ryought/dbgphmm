//!
//! Test of posterior sampling
//!
use super::super::{CopyNums, MultiDbg};
use super::Posterior;
use crate::e2e::Dataset;
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;
use crate::utils::timer;

///
/// 1. Create draft MultiDbg
/// 2. Infer posterior distribution
/// 3. Check the posterior probability of true copynums
/// 4. Check the posterior probability of copy num of edge of false-kmer
///
fn test_posterior(dataset: &Dataset, k: usize, param_infer: PHMMParams, gfa_filename: &str) {
    let (mut mdbg, t) = timer(|| MultiDbg::create_draft_from_dataset(k, &dataset));
    let paths_true = mdbg
        .compact_paths_from_styled_seqs(dataset.genome())
        .unwrap();
    // mdbg.to_paths_file("simple.paths", &paths_true);

    let post = mdbg.sample_posterior_with_dataset(&dataset, param_infer, 200, 10, 1, 10);
    // post.to_file("simple.post");
    let copy_nums_true = mdbg.copy_nums_from_compact_path(&paths_true);
    mdbg.set_copy_nums(&copy_nums_true);
    mdbg.to_gfa_post_file(gfa_filename, &post);

    println!("{}", mdbg.to_inspect_string(&post, &copy_nums_true));

    check_posterior_is_correct(&mdbg, &post, &copy_nums_true);
}

///
/// * P(true copy nums | R) is the highest
/// * P(X(e)=0 | R) is not too high for true edge e
///
fn check_posterior_is_correct(dbg: &MultiDbg, posterior: &Posterior, copy_nums_true: &CopyNums) {
    // (1)
    // for each edges whose true copy number is non-zero
    //
    for edge in dbg.graph_compact().edge_indices() {
        let copy_num_true = copy_nums_true[edge];
        if copy_num_true > 0 {
            assert!(posterior.p_edge_x(edge, 0) < Prob::from_prob(0.8));
        }
    }

    // (2)
    assert_eq!(posterior.max_copy_nums(), copy_nums_true);
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::e2e::{generate_dataset, ReadType};
    use crate::genome;

    #[test]
    #[ignore]
    fn simple_haploid() {
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

        test_posterior(&dataset, 20, PHMMParams::uniform(0.001), "shap.gfa");
    }

    #[test]
    #[ignore]
    fn simple_diploid() {
        let (genome, genome_size) = genome::diploid(1000, 0, 0.01, 0);
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

        test_posterior(&dataset, 20, PHMMParams::uniform(0.001), "sdip.gfa");
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

        test_posterior(&dataset, 20, PHMMParams::uniform(0.01), "repeat1k.gfa");
    }
}
