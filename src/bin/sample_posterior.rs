use clap::{AppSettings, ArgEnum, Parser};
use dbgphmm::common::collection::starts_and_ends_of_genome;
use dbgphmm::dbg::draft::EndNodeInference;
use dbgphmm::dbg::greedy::get_max_posterior_instance;
use dbgphmm::dbg::hashdbg_v2::HashDbg;
use dbgphmm::dbg::{Dbg, SimpleDbg};
use dbgphmm::e2e::{generate_dataset, Experiment, ReadType};
use dbgphmm::genome;
use dbgphmm::graph::cycle::CycleWithDir;
use dbgphmm::kmer::common::kmers_to_string_pretty;
use dbgphmm::kmer::VecKmer;
use dbgphmm::prelude::*;
use dbgphmm::utils::{check_memory_usage, timer};
use git_version::git_version;
use rayon::prelude::*;
use std::fs::File;
use std::io::prelude::*;
use std::path::PathBuf;

const GIT_VERSION: &str = git_version!();

#[derive(Parser, Debug)]
struct Opts {
    // dbg settings
    #[clap(long)]
    k_init: usize,
    #[clap(long)]
    k_final: usize,
    // read related
    #[clap(short = 'c')]
    coverage: usize,
    #[clap(short = 'p', default_value = "0.001")]
    p_error: f64,
    #[clap(short = 'l')]
    read_length: usize,
    #[clap(long)]
    use_fragment_read: bool,
    // genome related
    #[clap(short = 'U')]
    unit_size: usize,
    #[clap(short = 'N')]
    n_unit: usize,
    #[clap(short = 'E')]
    end_length: usize,
    #[clap(short = 'D', default_value = "0.01")]
    unit_divergence: f64,
    #[clap(short = 'H', default_value = "0.01")]
    hap_divergence: f64,
    #[clap(short = 'P', default_value = "1")]
    n_haplotypes: usize,
    // search related
    #[clap(long = "sigma", default_value = "100")]
    sigma: usize,
    #[clap(short = 's', default_value = "0")]
    seed: u64,
    #[clap(short = 'd', default_value = "10")]
    neighbor_depth: usize,
    #[clap(short = 'm', default_value = "3")]
    max_move: usize,
    /// purge threshold
    #[clap(long = "p0", default_value = "0.8")]
    p_0: f64,
    #[clap(long)]
    p_infer: Option<f64>,
    #[clap(long)]
    start_from_true: bool,
    /// use true dbg (k=k_init) as initial dbg and extend it until k reaches k_final
    #[clap(long)]
    start_from_true_dbg: bool,
    #[clap(long)]
    dbgviz_output: Option<PathBuf>,
    #[clap(long)]
    use_true_end_nodes: bool,
    /// construct true dbg (for each k) and infer the copy numbers on it.
    #[clap(long)]
    use_true_dbg: bool,
    #[clap(long)]
    use_homo_ends: bool,
    #[clap(long, default_value = "1")]
    copy_num_multiplicity: usize,
}

fn main() {
    let opts: Opts = Opts::parse();

    println!("# started_at={}", chrono::Local::now());
    println!("# version={}", GIT_VERSION);
    println!("# opts={:?}", opts);
    check_memory_usage();

    let (genome, genome_size) = if opts.use_homo_ends {
        genome::tandem_repeat_polyploid_with_unique_homo_ends(
            opts.unit_size,
            opts.n_unit,
            opts.seed,
            opts.end_length,
            opts.n_haplotypes,
            opts.hap_divergence,
            opts.seed,
        )
    } else {
        // will deprecate
        genome::tandem_repeat_polyploid_with_unique_ends(
            opts.unit_size,
            opts.n_unit,
            opts.unit_divergence,
            opts.seed,
            opts.seed,
            opts.end_length,
            opts.n_haplotypes,
            opts.hap_divergence,
            opts.seed,
        )
    };

    // let (genome, genome_size) = genome::tandem_repeat_diploid_example_ins();
    let coverage = opts.coverage;
    let param = PHMMParams::uniform(opts.p_error);
    let (dataset, mut dbg) = if opts.use_fragment_read {
        let dataset = generate_dataset(
            genome.clone(),
            genome_size,
            opts.seed, // read seed
            coverage,
            opts.read_length,
            ReadType::FragmentWithRevComp,
            param,
        );
        let end_node = if opts.use_true_end_nodes {
            EndNodeInference::Custom(starts_and_ends_of_genome(&genome, opts.k_init))
        } else {
            EndNodeInference::Auto
        };
        let dbg: SimpleDbg<VecKmer> =
            SimpleDbg::create_draft_from_fragment_seqs_with_adjusted_coverage(
                opts.k_init,
                dataset.reads(),
                dataset.coverage(),
                dataset.reads().average_length(),
                dataset.params().p_error().to_value(),
                &end_node,
            );
        (dataset, dbg)
    } else {
        let dataset = generate_dataset(
            genome.clone(),
            genome_size,
            opts.seed, // read seed
            coverage,
            genome_size * 2,
            ReadType::FullLength,
            param,
        );
        let dbg =
            SimpleDbg::create_draft_from_seqs(opts.k_init, dataset.reads(), dataset.coverage());
        (dataset, dbg)
    };

    if opts.start_from_true_dbg {
        dbg = SimpleDbg::from_styled_seqs(opts.k_init, dataset.genome());
    }

    // dataset.show_genome();
    // dataset.show_reads();
    dataset.show_reads_with_genome();
    let mut k = dbg.k();

    // phmm param used in inference and posterior-sampling.
    let param_infer = if let Some(p_infer) = opts.p_infer {
        PHMMParams::uniform(p_infer)
    } else {
        PHMMParams::uniform(opts.p_error)
    };
    check_memory_usage();

    while k <= opts.k_final {
        if opts.use_true_dbg {
            dbg = SimpleDbg::from_styled_seqs(k, dataset.genome());
        }
        let (copy_nums_true, _) = dbg
            .to_copy_nums_of_styled_seqs(&genome)
            .unwrap_or_else(|err| panic!("{}", err));
        if opts.start_from_true {
            let c = copy_nums_true.clone() * opts.copy_num_multiplicity;
            dbg.set_node_copy_nums(&c);
        }
        println!("# k={}", dbg.k());
        assert_eq!(dbg.k(), k);
        println!("# k={} n_dead_nodes={}", k, dbg.n_dead_nodes());
        println!("# k={} n_nodes={}", k, dbg.n_nodes());
        println!("# k={} n_edges={}", k, dbg.n_edges());
        println!("# k={} copy_num_stats={:?}", k, dbg.copy_num_stats());
        println!("# k={} degree_stats={:?}", k, dbg.degree_stats());
        let edbg = dbg.to_compact_edbg_graph();
        println!("# k={} n_nodes_compacted_edbg={}", k, edbg.node_count());
        println!("# k={} n_edges_compacted_edbg={}", k, edbg.edge_count());
        println!(
            "# k={} init_dist_from_true={}",
            k,
            dbg.to_node_copy_nums().dist(&copy_nums_true)
        );

        check_memory_usage();
        let (neighbors, t) =
            timer(|| dbg.neighbor_copy_nums_fast_compact_with_info(opts.neighbor_depth, false));
        eprintln!("[neighbors] N={} t={}", neighbors.len(), t);

        let (distribution, t) = timer(|| {
            dbg.search_posterior(
                dataset.reads(),
                param_infer,
                opts.neighbor_depth,
                opts.max_move,
                dataset.genome_size(),
                opts.sigma,
                |instance| {
                    println!("G\t{}\t{}", instance.info_string(), instance.move_count());
                },
            )
        });
        eprintln!("[search] k={} {}", k, t);

        println!(
            "#N\tk\tP(G|R)\tP(R|G)\tP(G)\tG\tn_haps\tmove_count\tdist_from_true\tmax_abs_diff_from_true\tcount_missing_and_error_kmers\tcycle_summary\tmissings\terrors\tdbg"
        );
        for (p_gr, instance, score) in distribution.iter() {
            dbg.set_node_copy_nums(instance.copy_nums());
            let ((n_missing, n_missing_null), (n_error, n_error_null)) = dbg.inspect_kmers(&genome);
            let (missings, errors) = dbg.missing_error_kmers(&genome);
            println!(
                "N\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                k,
                p_gr,
                score.p_rg().to_log_value(),
                score.p_g().to_log_value(),
                dbg.genome_size(),
                dbg.n_starting_kmers(),
                instance.move_count(),
                instance.copy_nums().dist(&copy_nums_true),
                instance.copy_nums().max_abs_diff(&copy_nums_true),
                format!(
                    "({:<3}{:<3}),({:<3}{:<3})",
                    n_missing, n_missing_null, n_error, n_error_null,
                ),
                instance.info_string(),
                kmers_to_string_pretty(&missings),
                kmers_to_string_pretty(&errors),
                dbg,
            );
        }

        // compare dense score of dbg_true and dbg_max
        // (a) dbg_true
        // dbg.set_node_copy_nums(&copy_nums_true);
        // let r_true = dbg.evaluate(param_infer, dataset.reads(), genome_size, opts.sigma);
        // println!("NT\t{}\t{}\t", k, r_true);
        // dbg.show_mapping_summary_for_reads(dataset.params(), dataset.reads());

        // (b) dbg_max
        // dbg.set_node_copy_nums(get_max_posterior_instance(&distribution).copy_nums());
        // let r_max = dbg.evaluate(param_infer, dataset.reads(), genome_size, opts.sigma);
        // println!("NM\t{}\t{}\t", k, r_max);
        // dbg.show_mapping_summary_for_reads(dataset.params(), dataset.reads());

        // set to max instance copy_nums in distribution
        dbg.set_node_copy_nums(get_max_posterior_instance(&distribution).copy_nums());
        let neighbors: Vec<_> = distribution
            .iter()
            .map(|(p_gr, instance, _score)| (instance.copy_nums().clone(), *p_gr))
            .collect();
        let read_count = HashDbg::from_seqs(k, dataset.reads());
        dbg.inspect_kmer_variance(&neighbors, &copy_nums_true, &read_count);

        if let Some(path) = &opts.dbgviz_output {
            let mut dbg_true = dbg.clone();
            dbg_true.set_node_copy_nums(&copy_nums_true);
            let dist = dbg.to_kmer_distribution(&neighbors);
            let copy_num_expected = dbg_true.to_copy_num_expected_vector(&dist);
            let json = dbg_true.to_cytoscape_with_info(
                |node| Some(format!("{}", dist[node.index()])),
                Some(&copy_num_expected),
            );
            let mut file = File::create(path.with_extension(format!("k{}.json", k))).unwrap();
            writeln!(file, "{}", json).unwrap();
        }

        let n_purged = dbg.purge_zero_copy_with_high_prob_kmer(
            &dbg.to_kmer_distribution(&neighbors),
            Prob::from_prob(opts.p_0),
        );
        println!("# k={} n_purged={}", dbg.k(), n_purged);

        // upgrade
        if opts.use_true_dbg {
            k += 1;
        } else {
            dbg = dbg.to_k_max_dbg_naive(opts.k_final);
            if k == dbg.k() {
                break;
            }
            k = dbg.k();
        }
    }

    println!("# finished_at={}", chrono::Local::now());
}
