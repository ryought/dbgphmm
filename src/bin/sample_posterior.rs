use clap::{AppSettings, ArgEnum, Clap};
use dbgphmm::e2e::{generate_experiment_with_draft, Experiment, ReadType};
use dbgphmm::genome;
use dbgphmm::graph::cycle::CycleWithDir;
use dbgphmm::kmer::common::kmers_to_string;
use dbgphmm::prelude::*;
use rayon::prelude::*;
use std::time::{Duration, Instant};

#[derive(Clap, Debug)]
struct Opts {
    #[clap(short = 'k', default_value = "40")]
    k: usize,
    #[clap(short = 'c', default_value = "20")]
    coverage: usize,
    #[clap(short = 'p', default_value = "0.001")]
    p_error: f64,
    #[clap(short = 'U', default_value = "100")]
    unit_size: usize,
    #[clap(short = 'N', default_value = "10")]
    n_unit: usize,
    #[clap(long = "div", default_value = "0.01")]
    unit_divergence: f64,
    #[clap(long = "sigma", default_value = "100")]
    sigma: usize,
    #[clap(short = 's', default_value = "0")]
    seed: u64,
    #[clap(short = 'd', default_value = "10")]
    neighbor_depth: usize,
    #[clap(short = 'm', default_value = "1")]
    max_move: usize,
}

fn main() {
    let opts: Opts = Opts::parse();
    let (genome, genome_size) = genome::tandem_repeat_haploid(
        opts.unit_size,
        opts.n_unit,
        opts.unit_divergence,
        opts.seed,
        opts.seed,
    );
    let coverage = opts.coverage;
    let param = PHMMParams::uniform(opts.p_error);
    let experiment = generate_experiment_with_draft(
        genome.clone(),
        genome_size,
        0,
        param,
        coverage,
        genome_size * 2,
        ReadType::FullLength,
        opts.k,
        40,
    );
    let mut dbg = experiment.dbg_draft_true.clone().unwrap();
    let (copy_nums_true, _) = dbg.to_copy_nums_of_styled_seqs(&genome).unwrap();

    println!("# opts={:?}", opts);
    for i in 0..genome.len() {
        println!("# genome[{}]={}", i, genome[i]);
    }
    println!("# k={}", dbg.k());
    println!("# n_dead_nodes={}", dbg.n_dead_nodes());
    println!("# n_nodes={}", dbg.n_nodes());
    println!("# n_edges={}", dbg.n_edges());
    println!("# copy_num_stats={:?}", dbg.copy_num_stats());
    println!("# degree_stats={:?}", dbg.degree_stats());
    let edbg = dbg.to_compact_edbg_graph();
    println!("# n_nodes_compacted_edbg={}", edbg.node_count());
    println!("# n_edges_compacted_edbg={}", edbg.edge_count());

    let distribution = dbg.search_posterior(
        experiment.dataset(),
        opts.neighbor_depth,
        opts.max_move,
        experiment.genome_size(),
        opts.sigma,
    );

    println!(
        // "#N genome_size\tcopy_nums_diff\tp\tmissing_and_error_kmers\tcycle_summary\tdbg\tcopy_nums"
        "#N P(G|R)\tP(R|G)\tP(G)\tG\tmove_count\tdist_from_true\tmissing_and_error_kmers\tcycle_summary\tdbg\tcopy_nums"
    );
    for (p_gr, instance, score) in distribution.iter() {
        dbg.set_node_copy_nums(instance.copy_nums());
        let ((n_missing, n_missing_null), (n_error, n_error_null)) = dbg.inspect_kmers(&genome);
        println!(
            "N {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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

    dbg.set_node_copy_nums(&copy_nums_true);
    let neighbors: Vec<_> = distribution
        .iter()
        .map(|(p_gr, instance, _score)| (instance.copy_nums().clone(), *p_gr))
        .collect();
    dbg.inspect_kmer_variance(&neighbors);
}
