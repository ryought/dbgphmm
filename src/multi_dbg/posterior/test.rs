//!
//! Test of posterior sampling
//!
use super::super::{neighbors::NeighborConfig, CopyNums, MultiDbg, Path};
use super::{infer_posterior_by_extension, Posterior};
use crate::common::{CopyNum, ReadCollection, Seq};
use crate::e2e::Dataset;
use crate::hmmv2::hint::Mappings;
use crate::hmmv2::params::PHMMParams;
use crate::prob::{p, Prob};
use crate::utils::timer;

///
/// 1. Create draft MultiDbg
/// 2. Infer posterior distribution
/// 3. Check the posterior probability of true copynums
/// 4. Check the posterior probability of copy num of edge of false-kmer
///
pub fn test_posterior(
    dataset: &Dataset,
    k: usize,
    param_infer: PHMMParams,
    gfa_filename: &str,
) -> (MultiDbg, Posterior, CopyNums) {
    println!("# started_at={}", chrono::Local::now());

    let (mut mdbg, _t) = timer(|| MultiDbg::create_draft_from_dataset(k, &dataset));
    let paths_true = mdbg.paths_from_styled_seqs(dataset.genome()).unwrap();
    // mdbg.to_paths_file("simple.paths", &paths_true);

    let mappings = mdbg.generate_mappings(param_infer, dataset.reads(), None);

    let (post, t) = timer(|| {
        mdbg.sample_posterior(
            param_infer,
            dataset.reads(),
            &mappings,
            dataset.genome_size(),
            200,
            NeighborConfig {
                max_cycle_size: 10,
                max_flip: 2,
                use_long_cycles: true,
                ignore_cycles_passing_terminal: true,
                use_reducers: true,
            },
            10,
            false,
        )
    });
    println!("sampled in {}ms", t);

    // post.to_file("simple.post");
    let copy_nums_true = mdbg.copy_nums_from_full_path(&paths_true);
    mdbg.set_copy_nums(&copy_nums_true);
    mdbg.to_gfa_post_file(gfa_filename, &post, Some(&copy_nums_true))
        .unwrap();
    mdbg.to_inspect_file(
        format!("{}.inspect", gfa_filename),
        &post,
        Some(&copy_nums_true),
    )
    .unwrap();

    println!("# finished_at={}", chrono::Local::now());

    (mdbg, post, copy_nums_true)
}

///
///
///
pub fn test_posterior_from_true<P: AsRef<std::path::Path>>(
    dataset: &Dataset,
    mut dbg: MultiDbg,
    sigma: CopyNum,
    inspect_filename: P,
    param: PHMMParams,
) -> (MultiDbg, Posterior, CopyNums) {
    let paths_true = dbg.paths_from_styled_seqs(dataset.genome()).unwrap();
    let copy_nums_true = dbg.copy_nums_from_full_path(&paths_true);
    dbg.set_copy_nums(&copy_nums_true);
    let mappings = dbg.generate_mappings(dataset.params(), dataset.reads(), None);

    let post = dbg.sample_posterior(
        param,
        dataset.reads(),
        &mappings,
        dataset.genome_size(),
        sigma,
        NeighborConfig {
            max_cycle_size: 10,
            // full=15 short=10
            max_flip: 2,
            use_long_cycles: true,
            ignore_cycles_passing_terminal: false,
            use_reducers: true,
        },
        10,
        false,
    );
    dbg.to_inspect_file(inspect_filename, &post, Some(&copy_nums_true))
        .unwrap();

    (dbg, post, copy_nums_true)
}

///
///
///
pub fn test_inference<P: AsRef<std::path::Path>>(
    dataset: &Dataset,
    k_init: usize,
    k_final: usize,
    param_infer: PHMMParams,
    param_error: PHMMParams,
    sigma: CopyNum,        // 200
    max_iter: usize,       // 10
    max_cycle_size: usize, // 10
    output_prefix: P,
) -> (MultiDbg, Posterior, Option<Vec<Path>>, Mappings) {
    let (dbg, _t) = timer(|| MultiDbg::create_draft_from_dataset(k_init, &dataset));
    test_inference_from_dbg_with_dataset(
        dataset,
        dbg,
        k_final,
        param_infer,
        param_error,
        sigma,
        max_iter,
        max_cycle_size,
        p(0.8),
        output_prefix,
        None,
    )
}

///
/// Mapping extension test
///
/// Compare
/// * Exact mapping calculated by forward/backward sparse
/// * Inherited mapping by extension
///
pub fn test_mapping_extension<P: AsRef<std::path::Path>>(
    dataset: &Dataset,
    dbg: MultiDbg,
    k_final: usize,
    param_infer: PHMMParams, // virtual error rate
    output_prefix: P,
) {
    let mut dbg = dbg;
    let mut paths = dbg.paths_from_styled_seqs(dataset.genome()).ok();
    let mut mappings = dbg.generate_mappings(param_infer, dataset.reads(), None);
    let output: std::path::PathBuf = output_prefix.as_ref().into();

    eprintln!("generating mappings");

    while dbg.k() < k_final {
        eprintln!("k={}", dbg.k());

        // compute mapping by extension and refine
        eprintln!("computing extend...");
        let copy_nums = dbg.copy_nums_from_full_path(&paths.as_ref().unwrap());
        eprintln!("copy_nums={}", copy_nums);
        dbg.set_copy_nums(&copy_nums);
        let zero_edges: Vec<_> = dbg
            .graph_compact()
            .edge_indices()
            .filter(|&e| dbg.copy_num_of_edge_in_compact(e) == 0)
            .collect();
        eprintln!("zero_edges={:?}", zero_edges);
        let t_start_extend = std::time::Instant::now();
        eprintln!("extending..");
        eprintln!("genome_size={}", dbg.genome_size());
        let (dbg_new, paths_new, mappings_new) =
            dbg.purge_and_extend(&zero_edges, k_final, true, paths, Some(&mappings));
        dbg = dbg_new;
        paths = paths_new;
        mappings = mappings_new.unwrap();
        eprintln!("genome_size2={}", dbg.genome_size());
        let t_extend = t_start_extend.elapsed();
        eprintln!("extend t={}ms", t_extend.as_millis());
        let t_start_refine = std::time::Instant::now();
        eprintln!("refining..");
        mappings = dbg.generate_mappings(param_infer, dataset.reads(), Some(&mappings));
        let t_refine = t_start_refine.elapsed();
        eprintln!("refine t={}ms", t_refine.as_millis());

        dbg.to_map_file(
            output.with_extension(format!("k{}.extend.map", dbg.k())),
            dataset.reads(),
            &mappings,
        )
        .unwrap();

        // compute true mapping
        eprintln!("computing true...");
        let t_start_map = std::time::Instant::now();
        let mappings_true = dbg.generate_mappings(param_infer, dataset.reads(), None);
        let t_map = t_start_map.elapsed();
        eprintln!("map t={}ms", t_map.as_millis());

        dbg.to_map_file(
            output.with_extension(format!("k{}.true.map", dbg.k())),
            dataset.reads(),
            &mappings_true,
        )
        .unwrap();

        dbg.to_dbg_file(output.with_extension(format!("k{}.dbg", dbg.k())))
            .unwrap();
        dbg.to_gfa_file(output.with_extension(format!("k{}.gfa", dbg.k())))
            .unwrap();

        let p_extend = dbg.to_likelihood(param_infer, dataset.reads(), Some(&mappings));
        let p_true = dbg.to_likelihood(param_infer, dataset.reads(), Some(&mappings_true));
        let p_true2 = dbg.to_likelihood(param_infer, dataset.reads(), None);

        let phmm = dbg.to_uniform_phmm(param_infer);
        let p_unif_true = phmm.to_full_prob_reads(dataset.reads(), Some(&mappings_true), true);
        let p_unif_true2 = phmm.to_full_prob_reads(dataset.reads(), None, true);

        println!(
            "k={} p_extend={} p_true={} p_true2={} p_unif_true={} p_unif_true2={}  t_extend={} t_refine={} t_map={}",
            dbg.k(),
            p_extend,
            p_true,
            p_true2,
            p_unif_true,
            p_unif_true2,
            t_extend.as_millis(),
            t_refine.as_millis(),
            t_map.as_millis()
        );
    }
}

///
///
///
pub fn test_inference_from_dbg<S: Seq, P: AsRef<std::path::Path>>(
    dbg: MultiDbg,
    reads: &ReadCollection<S>,
    k_final: usize,
    param_infer: PHMMParams, // virtual error rate
    param_error: PHMMParams, // actual error rate of reads
    genome_size: CopyNum,
    genome_size_sigma: CopyNum, // 200
    max_iter: usize,            // 10
    max_cycle_size: usize,      // 10
    p0: Prob,
    output_prefix: P,
    paths: Option<Vec<Path>>,
    mappings: Option<Mappings>,
) -> (MultiDbg, Posterior, Option<Vec<Path>>, Mappings) {
    println!("# started_at={}", chrono::Local::now());

    let output: String = output_prefix.as_ref().to_str().unwrap().to_owned();

    let (dbg, posterior, paths, mappings) = infer_posterior_by_extension(
        k_final,
        dbg,
        param_infer,
        param_error,
        reads,
        genome_size,
        genome_size_sigma,
        NeighborConfig {
            max_cycle_size,
            max_flip: 2,
            use_long_cycles: true,
            ignore_cycles_passing_terminal: true,
            use_reducers: true,
        },
        max_iter,
        p0,
        |dbg, posterior, paths, mappings| {
            let k = dbg.k();
            println!("callback k={} n_edges={}", k, dbg.n_edges_full());

            dbg.to_dbg_file(format!("{}.k{}.dbg", output, k)).unwrap();
            posterior
                .to_file(format!("{}.k{}.post", output, k))
                .unwrap();

            let copy_nums_true = paths
                .as_ref()
                .map(|paths| dbg.copy_nums_from_full_path(paths));
            dbg.to_gfa_post_file(
                format!("{}.k{}.gfa", output, k),
                posterior,
                copy_nums_true.as_ref(),
            )
            .unwrap();
            dbg.to_inspect_file(
                format!("{}.k{}.inspect", output, k),
                posterior,
                copy_nums_true.as_ref(),
            )
            .unwrap();
            dbg.to_map_file(format!("{}.k{}.mpz", output, k), reads, mappings)
                .unwrap();
        },
        |_dbg, _mappings| {},
        |_dbg, _mappings| {},
        paths,
        mappings,
    );

    // output final
    let copy_nums_true = paths
        .as_ref()
        .map(|paths| dbg.copy_nums_from_full_path(paths));
    dbg.to_dbg_file(format!("{}.final.dbg", output)).unwrap();
    dbg.to_gfa_post_file(
        format!("{}.final.gfa", output),
        &posterior,
        copy_nums_true.as_ref(),
    )
    .unwrap();
    dbg.to_inspect_file(
        format!("{}.final.inspect", output),
        &posterior,
        copy_nums_true.as_ref(),
    )
    .unwrap();
    dbg.to_fasta_linear(format!("{}.final.euler.fa", output))
        .unwrap();

    println!("# finished_at={}", chrono::Local::now());

    (dbg, posterior, paths, mappings)
}

///
///
///
pub fn test_inference_from_dbg_with_dataset<P: AsRef<std::path::Path>>(
    dataset: &Dataset,
    dbg: MultiDbg,
    k_final: usize,
    param_infer: PHMMParams, // virtual error rate
    param_error: PHMMParams, // actual error rate of reads
    sigma: CopyNum,          // 200
    max_iter: usize,         // 10
    max_cycle_size: usize,   // 10
    p0: Prob,
    output_prefix: P,
    mappings: Option<Mappings>,
) -> (MultiDbg, Posterior, Option<Vec<Path>>, Mappings) {
    let paths_true = dbg.paths_from_styled_seqs(dataset.genome());
    test_inference_from_dbg(
        dbg,
        dataset.reads(),
        k_final,
        param_infer,
        param_error,
        dataset.genome_size(),
        sigma,
        max_iter,
        max_cycle_size,
        p0,
        output_prefix,
        paths_true.ok(),
        mappings,
    )
}

///
/// * P(X(e)=0 | R) is not too high for true edge e
///
#[allow(dead_code)]
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
#[allow(dead_code)]
fn check_posterior_highest_at_true(
    _dbg: &MultiDbg,
    posterior: &Posterior,
    copy_nums_true: &CopyNums,
) {
    // (2)
    assert_eq!(posterior.max_copy_nums(), copy_nums_true);
}

#[allow(dead_code)]
fn check_posterior_highest_at_true_path(
    dbg: &MultiDbg,
    posterior: &Posterior,
    paths_true: &Option<Vec<Path>>,
) {
    assert!(paths_true.is_some());
    let copy_nums_true = dbg.copy_nums_from_full_path(paths_true.as_ref().unwrap());
    assert_eq!(posterior.max_copy_nums(), &copy_nums_true);
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
        let genome = genome::simple(1000, 0);
        let param = PHMMParams::uniform(0.005);
        let dataset = generate_dataset(genome, 0, 20, 500, ReadType::FragmentWithRevComp, param);
        dataset.show_reads_with_genome();

        // without hints
        let (mdbg, post, copy_nums_true) =
            test_posterior(&dataset, 20, PHMMParams::uniform(0.001), "shap.gfa");
        check_posterior_non_zero_edges(&mdbg, &post, &copy_nums_true);
        check_posterior_highest_at_true(&mdbg, &post, &copy_nums_true);
    }

    #[test]
    #[ignore]
    fn simple_diploid_posterior() {
        let genome = genome::diploid(1000, 0, 0.01, 0);
        let param = PHMMParams::uniform(0.005);
        let dataset = generate_dataset(genome, 0, 20, 500, ReadType::FragmentWithRevComp, param);
        dataset.show_reads_with_genome();

        // if without hint, 35sec
        // if with hint, 10sec
        let (mdbg, post, copy_nums_true) =
            test_posterior(&dataset, 20, PHMMParams::uniform(0.001), "sdip.gfa");
        check_posterior_non_zero_edges(&mdbg, &post, &copy_nums_true);
        check_posterior_highest_at_true(&mdbg, &post, &copy_nums_true);
    }

    #[test]
    #[ignore]
    fn simple_diploid_inference() {
        let genome = genome::diploid(1000, 0, 0.01, 0);
        let param = PHMMParams::uniform(0.005);
        let dataset = generate_dataset(genome, 0, 20, 500, ReadType::FragmentWithRevComp, param);
        dataset.show_reads_with_genome();

        let (dbg, post, paths, _) = test_inference(
            &dataset,
            20,
            500,
            PHMMParams::uniform(0.001),
            PHMMParams::uniform(0.005),
            200,
            10,
            10,
            "sdip/p0001",
        );
        check_posterior_highest_at_true_path(&dbg, &post, &paths);
    }

    #[test]
    #[ignore]
    fn repeat_1k_posterior() {
        let genome = genome::tandem_repeat_polyploid_with_unique_homo_ends(
            100, 10, 0, 0.0, 0, 300, 2, 0.01, 0,
        );
        let param = PHMMParams::uniform(0.01);
        let dataset = generate_dataset(genome, 0, 20, 500, ReadType::FragmentWithRevComp, param);
        dataset.show_reads_with_genome();

        // if without hint, 35sec
        // if with hint, 10sec
        let (mdbg, post, copy_nums_true) =
            test_posterior(&dataset, 20, PHMMParams::uniform(0.001), "repeat1k.gfa");
        check_posterior_non_zero_edges(&mdbg, &post, &copy_nums_true);
        // below is not satisfied in tandem repeat
        // check_posterior_highest_at_true(&mdbg, &post, &copy_nums_true);
    }

    // #[test]
    #[ignore]
    fn repeat_1k_inference() {
        let genome = genome::tandem_repeat_polyploid_with_unique_homo_ends(
            100, 10, 0, 0.0, 0, 300, 2, 0.01, 0,
        );
        let param = PHMMParams::uniform(0.01);
        let dataset = generate_dataset(genome, 0, 20, 500, ReadType::FragmentWithRevComp, param);
        dataset.show_reads_with_genome();

        test_inference(
            &dataset,
            20,
            500,
            PHMMParams::uniform(0.001),
            PHMMParams::uniform(0.005),
            200,
            10,
            10,
            "repeat1k/p0001",
        );
    }

    // #[test]
    #[ignore]
    fn repeat_u200_inference() {
        let genome = genome::tandem_repeat_polyploid_with_unique_homo_ends(
            200, 10, 0, 0.0, 0, 200, 2, 0.01, 0,
        );
        let param = PHMMParams::uniform(0.01);
        let dataset = generate_dataset(genome, 0, 20, 500, ReadType::FragmentWithRevComp, param);
        dataset.show_reads_with_genome();

        test_inference(
            &dataset,
            20,
            500,
            PHMMParams::uniform(0.001),
            PHMMParams::uniform(0.005),
            200,
            10,
            10,
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
        let mappings = mdbg.generate_mappings(param, &reads, None);
        let hint = &mappings[0];
        println!("hint={:?}", hint);
        // most probable node for bases[i] is v[i+1]
        assert_eq!(hint.nodes(0)[0], ni(1));
        assert_eq!(hint.nodes(1)[0], ni(2));
        assert_eq!(hint.nodes(2)[0], ni(3));
        assert_eq!(hint.nodes(3)[0], ni(4));
        assert_eq!(hint.nodes(4)[0], ni(5));

        // (2) read from repetitive
        let reads = ReadCollection::from(vec![b"GCAGCAGG".to_vec()]);
        let mappings = mdbg.generate_mappings(param, &reads, None);
        let hint = &mappings[0];
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
