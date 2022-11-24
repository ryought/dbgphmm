use clap::{AppSettings, ArgEnum, Clap};
use dbgphmm::dbg::greedy::get_max_posterior_instance;
use dbgphmm::dbg::{Dbg, SimpleDbg};
use dbgphmm::e2e::{generate_dataset, Experiment, ReadType};
use dbgphmm::genome;
use dbgphmm::graph::cycle::CycleWithDir;
use dbgphmm::kmer::common::kmers_to_string;
use dbgphmm::kmer::VecKmer;
use dbgphmm::prelude::*;
use git_version::git_version;
use rayon::prelude::*;
use std::time::{Duration, Instant};

const GIT_VERSION: &str = git_version!();

#[derive(Clap, Debug)]
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
    #[clap(long)]
    start_from_true: bool,
}

fn main() {
    let opts: Opts = Opts::parse();

    println!("# started_at={}", chrono::Local::now());
    println!("# version={}", GIT_VERSION);
    println!("# opts={:?}", opts);

    let (genome, genome_size) = genome::tandem_repeat_polyploid_with_unique_ends(
        opts.unit_size,
        opts.n_unit,
        opts.unit_divergence,
        opts.seed,
        opts.seed,
        opts.end_length,
        opts.n_haplotypes,
        opts.hap_divergence,
        opts.seed,
    );
    let coverage = opts.coverage;
    let param = PHMMParams::uniform(opts.p_error);
    let (dataset, mut dbg) = if opts.use_fragment_read {
        let dataset = generate_dataset(
            genome.clone(),
            genome_size,
            0, // read seed
            coverage,
            opts.read_length,
            ReadType::FragmentWithRevComp,
            param,
        );
        let dbg: SimpleDbg<VecKmer> =
            SimpleDbg::create_draft_from_fragment_seqs_with_adjusted_coverage(
                opts.k_init,
                dataset.reads(),
                dataset.coverage(),
                dataset.reads().average_length(),
                dataset.params().p_error().to_value(),
            );
        (dataset, dbg)
    } else {
        let dataset = generate_dataset(
            genome.clone(),
            genome_size,
            0, // read seed
            coverage,
            genome_size * 2,
            ReadType::FullLength,
            param,
        );
        let dbg =
            SimpleDbg::create_draft_from_seqs(opts.k_init, dataset.reads(), dataset.coverage());
        (dataset, dbg)
    };

    if opts.start_from_true {
        let (copy_nums_true, _) = dbg
            .to_copy_nums_of_styled_seqs(&genome)
            .unwrap_or_else(|err| panic!("{}", err));
        dbg.set_node_copy_nums(&copy_nums_true);
    }

    dataset.show_genome();
    dataset.show_reads();
    let mut k = dbg.k();

    while k <= opts.k_final {
        let (copy_nums_true, _) = dbg
            .to_copy_nums_of_styled_seqs(&genome)
            .unwrap_or_else(|err| panic!("{}", err));
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

        let distribution = dbg.search_posterior(
            &dataset,
            opts.neighbor_depth,
            opts.max_move,
            dataset.genome_size(),
            opts.sigma,
        );

        println!(
            "#N\tk\tP(G|R)\tP(R|G)\tP(G)\tG\tmove_count\tdist_from_true\tmissing_and_error_kmers\tcycle_summary\tdbg\tcopy_nums"
        );
        for (p_gr, instance, score) in distribution.iter() {
            dbg.set_node_copy_nums(instance.copy_nums());
            let ((n_missing, n_missing_null), (n_error, n_error_null)) = dbg.inspect_kmers(&genome);
            println!(
                "N\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                k,
                p_gr,
                score.p_rg,
                score.p_g,
                dbg.genome_size(),
                instance.move_count(),
                instance.copy_nums().dist(&copy_nums_true),
                format!(
                    "({:<3}{:<3}),({:<3}{:<3})",
                    n_missing, n_missing_null, n_error, n_error_null,
                ),
                instance.info(),
                dbg,
                instance.copy_nums(),
            );
        }

        // set to max instance copy_nums in distribution
        dbg.set_node_copy_nums(get_max_posterior_instance(&distribution).copy_nums());
        let neighbors: Vec<_> = distribution
            .iter()
            .map(|(p_gr, instance, _score)| (instance.copy_nums().clone(), *p_gr))
            .collect();
        dbg.inspect_kmer_variance(&neighbors, &copy_nums_true);
        let n_purged = dbg.purge_zero_copy_with_high_prob_kmer(
            &dbg.to_kmer_distribution(&neighbors),
            Prob::from_prob(0.8),
        );
        println!("# k={} n_purged={}", dbg.k(), n_purged);

        // upgrade
        dbg = dbg.to_k_max_dbg_naive(opts.k_final);
        if k == dbg.k() {
            break;
        }
        k = dbg.k();
    }

    println!("# finished_at={}", chrono::Local::now());
}