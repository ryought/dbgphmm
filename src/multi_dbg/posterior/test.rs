//!
//! Test of posterior sampling
//!
use super::super::{CopyNums, MultiDbg};
use super::{infer_posterior_by_extension, Posterior};
use crate::e2e::Dataset;
use crate::hmmv2::params::PHMMParams;
use crate::prob::{p, Prob};
use crate::utils::timer;

///
/// 1. Create draft MultiDbg
/// 2. Infer posterior distribution
/// 3. Check the posterior probability of true copynums
/// 4. Check the posterior probability of copy num of edge of false-kmer
///
fn test_posterior(
    dataset: &Dataset,
    k: usize,
    param_infer: PHMMParams,
    gfa_filename: &str,
    use_hint: bool,
) -> (MultiDbg, Posterior, CopyNums) {
    println!("# started_at={}", chrono::Local::now());

    let (mut mdbg, t) = timer(|| MultiDbg::create_draft_from_dataset(k, &dataset));
    let paths_true = mdbg.paths_from_styled_seqs(dataset.genome()).unwrap();
    // mdbg.to_paths_file("simple.paths", &paths_true);

    let reads = if use_hint {
        mdbg.generate_hints(param_infer, dataset.reads().clone(), true)
    } else {
        dataset.reads().clone()
    };

    let (post, t) = timer(|| {
        mdbg.sample_posterior(
            param_infer,
            &reads,
            dataset.genome_size(),
            200,
            10,
            1,
            10,
            true,
        )
    });
    println!("sampled in {}ms", t);

    // post.to_file("simple.post");
    let copy_nums_true = mdbg.copy_nums_from_full_path(&paths_true);
    mdbg.set_copy_nums(&copy_nums_true);
    mdbg.to_gfa_post_file(gfa_filename, &post);
    mdbg.to_inspect_file(
        format!("{}.inspect", gfa_filename),
        &post,
        Some(&copy_nums_true),
    );

    println!("# finished_at={}", chrono::Local::now());

    (mdbg, post, copy_nums_true)
}

///
///
///
fn test_inference(
    dataset: &Dataset,
    k_init: usize,
    k_final: usize,
    param_infer: PHMMParams,
    output_prefix: &str,
) {
    println!("# started_at={}", chrono::Local::now());

    let (dbg, t) = timer(|| MultiDbg::create_draft_from_dataset(k_init, &dataset));
    let paths_true = dbg.paths_from_styled_seqs(dataset.genome());
    let reads = dbg.generate_hints(param_infer, dataset.reads().clone(), true);
    let output: std::path::PathBuf = output_prefix.into();

    infer_posterior_by_extension(
        k_final,
        dbg,
        param_infer,
        reads,
        dataset.genome_size(),
        200,
        10,
        1,
        10,
        p(0.8),
        |dbg, posterior, paths, reads| {
            let k = dbg.k();
            println!("callback k={} n_edges={}", k, dbg.n_edges_full());

            dbg.to_dbg_file(output.with_extension(format!("k{}.dbg", k)));
            posterior.to_file(output.with_extension(format!("k{}.post", k)));
            dbg.to_gfa_post_file(output.with_extension(format!("k{}.gfa", k)), posterior);

            let copy_nums_true = paths
                .as_ref()
                .map(|paths| dbg.copy_nums_from_full_path(paths));
            dbg.to_inspect_file(
                output.with_extension(format!("k{}.inspect", k)),
                posterior,
                copy_nums_true.as_ref(),
            );
        },
        paths_true.ok(),
    );

    println!("# finished_at={}", chrono::Local::now());
}

///
/// * P(X(e)=0 | R) is not too high for true edge e
///
fn check_posterior_non_zero_edges(
    dbg: &MultiDbg,
    posterior: &Posterior,
    copy_nums_true: &CopyNums,
) {
    // (1)
    // for each edges whose true copy number is non-zero
    //
    for edge in dbg.graph_compact().edge_indices() {
        let copy_num_true = copy_nums_true[edge];
        if copy_num_true > 0 {
            assert!(posterior.p_edge_x(edge, 0) < Prob::from_prob(0.8));
        }
    }
}

///
/// * P(true copy nums | R) is the highest
///
fn check_posterior_highest_at_true(
    dbg: &MultiDbg,
    posterior: &Posterior,
    copy_nums_true: &CopyNums,
) {
    // (2)
    assert_eq!(posterior.max_copy_nums(), copy_nums_true);
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::super::super::toy;
    use super::*;
    use crate::common::{ni, ReadCollection, Seq};
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

        // without hints
        let (mdbg, post, copy_nums_true) =
            test_posterior(&dataset, 20, PHMMParams::uniform(0.001), "shap.gfa", false);
        check_posterior_non_zero_edges(&mdbg, &post, &copy_nums_true);
        check_posterior_highest_at_true(&mdbg, &post, &copy_nums_true);

        // with hints
        let (mdbg, post, copy_nums_true) =
            test_posterior(&dataset, 20, PHMMParams::uniform(0.001), "shap.gfa", true);
        check_posterior_non_zero_edges(&mdbg, &post, &copy_nums_true);
        check_posterior_highest_at_true(&mdbg, &post, &copy_nums_true);
    }

    #[test]
    #[ignore]
    fn simple_diploid_posterior() {
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

        // if without hint, 35sec
        // if with hint, 10sec
        let (mdbg, post, copy_nums_true) =
            test_posterior(&dataset, 20, PHMMParams::uniform(0.001), "sdip.gfa", true);
        check_posterior_non_zero_edges(&mdbg, &post, &copy_nums_true);
        check_posterior_highest_at_true(&mdbg, &post, &copy_nums_true);
    }

    #[test]
    #[ignore]
    fn simple_diploid_inference() {
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

        test_inference(&dataset, 20, 100, PHMMParams::uniform(0.001), "sdip/sdip");
    }

    #[test]
    #[ignore]
    fn repeat_1k_posterior() {
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

        // if without hint, 35sec
        // if with hint, 10sec
        let (mdbg, post, copy_nums_true) = test_posterior(
            &dataset,
            20,
            PHMMParams::uniform(0.001),
            "repeat1k.gfa",
            true,
        );
        check_posterior_non_zero_edges(&mdbg, &post, &copy_nums_true);
        // below is not satisfied in tandem repeat
        // check_posterior_highest_at_true(&mdbg, &post, &copy_nums_true);
    }

    #[test]
    #[ignore]
    fn repeat_1k_inference() {
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

        test_inference(
            &dataset,
            20,
            500,
            PHMMParams::uniform(0.001),
            "repeat1k/p0001",
        );
    }

    #[test]
    #[ignore]
    fn repeat_u200_inference() {
        let (genome, genome_size) =
            genome::tandem_repeat_polyploid_with_unique_homo_ends(200, 10, 0, 200, 2, 0.01, 0);
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

        test_inference(
            &dataset,
            20,
            500,
            PHMMParams::uniform(0.001),
            "repeat_u200/p0001",
        );
    }

    #[test]
    fn hint_for_toy() {
        let mdbg = toy::repeat();
        mdbg.show_graph_with_kmer();
        let param = PHMMParams::uniform(0.01);
        let phmm = mdbg.to_uniform_phmm(param);
        println!("{}", phmm);

        // (1) read from first
        let reads = ReadCollection::from(vec![b"CCCAG".to_vec()]);
        let reads_with_hint = mdbg.generate_hints(param, reads, false);
        for (r, h) in reads_with_hint.iter_with_hint() {
            println!("{} {:?}", r.to_str(), h);
        }
        let hint = reads_with_hint.hint(0);
        // most probable node for bases[i] is v[i+1]
        assert_eq!(hint.nodes(0)[0], ni(1));
        assert_eq!(hint.nodes(1)[0], ni(2));
        assert_eq!(hint.nodes(2)[0], ni(3));
        assert_eq!(hint.nodes(3)[0], ni(4));
        assert_eq!(hint.nodes(4)[0], ni(5));

        // (2) read from repetitive
        let reads = ReadCollection::from(vec![b"GCAGCAGG".to_vec()]);
        let reads_with_hint = mdbg.generate_hints(param, reads, false);
        for (r, h) in reads_with_hint.iter_with_hint() {
            println!("{} {:?}", r.to_str(), h);
        }
        let hint = reads_with_hint.hint(0);
        assert_eq!(hint.nodes(0)[0], ni(8));
        assert_eq!(hint.nodes(1)[0], ni(6));
        assert_eq!(hint.nodes(2)[0], ni(7));
        assert_eq!(hint.nodes(3)[0], ni(8));
        assert_eq!(hint.nodes(4)[0], ni(6));
        assert_eq!(hint.nodes(5)[0], ni(7));
        assert_eq!(hint.nodes(6)[0], ni(8));
        assert_eq!(hint.nodes(7)[0], ni(9));
    }
}
