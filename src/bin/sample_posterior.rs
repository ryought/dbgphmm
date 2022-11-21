use clap::{AppSettings, ArgEnum, Clap};
use dbgphmm::dbg::greedy::get_max_posterior_instance;
use dbgphmm::e2e::{generate_experiment_with_draft, Experiment, ReadType};
use dbgphmm::genome;
use dbgphmm::graph::cycle::CycleWithDir;
use dbgphmm::kmer::common::kmers_to_string;
use dbgphmm::prelude::*;
use rayon::prelude::*;
use std::time::{Duration, Instant};

#[derive(Clap, Debug)]
struct Opts {
    #[clap(long, default_value = "40")]
    k_init: usize,
    #[clap(long, default_value = "50")]
    k_final: usize,
    #[clap(short = 'c', default_value = "20")]
    coverage: usize,
    #[clap(short = 'p', default_value = "0.001")]
    p_error: f64,
    #[clap(short = 'U', default_value = "100")]
    unit_size: usize,
    #[clap(short = 'N', default_value = "10")]
    n_unit: usize,
    #[clap(short = 'E', default_value = "50")]
    end_length: usize,
    #[clap(short = 'D', default_value = "0.01")]
    unit_divergence: f64,
    #[clap(short = 'H', default_value = "0.01")]
    hap_divergence: f64,
    #[clap(short = 'P', default_value = "1")]
    n_haplotypes: usize,
    #[clap(long = "sigma", default_value = "100")]
    sigma: usize,
    #[clap(short = 's', default_value = "0")]
    seed: u64,
    #[clap(short = 'd', default_value = "10")]
    neighbor_depth: usize,
    #[clap(short = 'm', default_value = "3")]
    max_move: usize,
    #[clap(long)]
    from_approx: bool,
}

fn main() {
    let opts: Opts = Opts::parse();
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
    let experiment = generate_experiment_with_draft(
        genome.clone(),
        genome_size,
        0, // read seed
        param,
        coverage,
        genome_size * 2,
        ReadType::FullLength,
        opts.k_init,
        40, // XXX this argument is to be removed
    );
    let mut dbg = if opts.from_approx {
        experiment.dbg_draft.clone().unwrap()
    } else {
        experiment.dbg_draft_true.clone().unwrap()
    };

    println!("# started_at={}", chrono::Local::now());
    println!("# opts={:?}", opts);
    println!("# n_hap={}", genome.len());
    for i in 0..genome.len() {
        println!("# genome[{}]={}", i, genome[i]);
    }
    let mut k = dbg.k();

    while k <= opts.k_final {
        let (copy_nums_true, _) = dbg.to_copy_nums_of_styled_seqs(&genome).unwrap();
        println!("# k={}", dbg.k());
        assert_eq!(dbg.k(), k);
        println!("# n_dead_nodes={}", dbg.n_dead_nodes());
        println!("# n_nodes={}", dbg.n_nodes());
        println!("# n_edges={}", dbg.n_edges());
        println!("# copy_num_stats={:?}", dbg.copy_num_stats());
        println!("# degree_stats={:?}", dbg.degree_stats());
        let edbg = dbg.to_compact_edbg_graph();
        println!("# n_nodes_compacted_edbg={}", edbg.node_count());
        println!("# n_edges_compacted_edbg={}", edbg.edge_count());
        println!(
            "# init_dist_from_true={}",
            dbg.to_node_copy_nums().dist(&copy_nums_true)
        );

        let distribution = dbg.search_posterior(
            experiment.dataset(),
            opts.neighbor_depth,
            opts.max_move,
            experiment.genome_size(),
            opts.sigma,
        );

        println!(
            "#N k\tP(G|R)\tP(R|G)\tP(G)\tG\tmove_count\tdist_from_true\tmissing_and_error_kmers\tcycle_summary\tdbg\tcopy_nums"
        );
        for (p_gr, instance, score) in distribution.iter() {
            dbg.set_node_copy_nums(instance.copy_nums());
            let ((n_missing, n_missing_null), (n_error, n_error_null)) = dbg.inspect_kmers(&genome);
            println!(
                "N {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
        dbg.inspect_kmer_variance(&neighbors);
        let n_purged = dbg.purge_zero_copy_with_high_prob_kmer(
            &dbg.to_kmer_distribution(&neighbors),
            Prob::from_prob(0.8),
        );
        println!("# k={} n_purged={}", dbg.k(), n_purged);

        // upgrade
        dbg = dbg.to_k_max_dbg_naive(opts.k_final);
        k = dbg.k();
    }

    println!("# finished_at={}", chrono::Local::now());
}
