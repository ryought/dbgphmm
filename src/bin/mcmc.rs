use clap::{AppSettings, ArgEnum, Clap};
use dbgphmm::e2e::{generate_experiment, Experiment, ReadType};
use dbgphmm::genome;
use dbgphmm::graph::cycle::CycleWithDir;
use dbgphmm::prelude::*;
use rayon::prelude::*;

#[derive(Clap)]
struct Opts {
    #[clap(short = 'k', default_value = "24")]
    k: usize,
    #[clap(short = 'c', default_value = "20")]
    coverage: usize,
    #[clap(short = 'p', default_value = "0.003")]
    p_error: f64,
    #[clap(short, long)]
    simple: bool,
}

fn run_mcmc() {
    let opts: Opts = Opts::parse();
    let (genome, genome_size) = if opts.simple {
        genome::simple(100, 5)
    } else {
        genome::tandem_repeat_haploid(20, 5, 0.01, 0, 0)
        // genome::tandem_repeat_haploid(50, 4, 0.05, 0, 0);
    };
    let coverage = opts.coverage;
    let param = PHMMParams::uniform(opts.p_error);
    let dataset = generate_experiment(
        genome.clone(),
        genome_size,
        0,
        param,
        coverage,
        2000,
        ReadType::FullLengthForHaploid,
        opts.k,
        40,
    );
    let dbg_raw = dataset.dbg_raw.clone();
    let mut dbg_true = dbg_raw.clone().shrink_single_copy_node();
    dbg_true.remove_zero_copy_node();
    let (copy_nums_true, _) = dbg_true.to_copy_nums_of_styled_seqs(&genome).unwrap();
    dbg_true.set_node_copy_nums(&copy_nums_true);

    println!("# genome={}", genome[0]);
    println!("# k={}", dbg_true.k());
    println!("# n_kmers_with_null={}", dbg_true.n_kmers_with_null());
    println!("# n_dead_nodes={}", dbg_true.n_dead_nodes());
    println!("# n_nodes={}", dbg_true.n_nodes());
    println!("# n_edges={}", dbg_true.n_edges());

    let mut neighbors = dbg_true.neighbor_copy_nums_and_cycles();
    println!("# n_neighbors={}", neighbors.len());
    neighbors.push((copy_nums_true.clone(), CycleWithDir::empty()));

    println!(
        "#N genome_size\tcopy_nums_diff\tp\tmissing_and_error_kmers\tcycle_summary\tdbg\tcopy_nums"
    );

    let mut neighbors: Vec<_> = neighbors
        .into_par_iter()
        .filter(|(copy_nums, _)| copy_nums.sum() > 0)
        .map(|(copy_nums, cycle)| {
            let mut dbg = dbg_true.clone();
            dbg.set_node_copy_nums(&copy_nums);
            let p = dbg.to_full_prob(param, dataset.reads());
            (copy_nums, cycle, dbg, p)
        })
        .collect();
    neighbors.sort_by_key(|(_, _, _, p)| *p);
    neighbors.reverse();
    neighbors.iter().for_each(|(copy_nums, cycle, dbg, p)| {
        let ((n_missing, n_missing_null), (n_error, n_error_null)) = dbg.inspect_kmers(&genome);
        println!(
            "N {:<100}{}",
            format!(
                "{:<5}\t{:<5}\t{:<20}\t{:<20}\t{:<30}",
                dbg.genome_size(),
                copy_nums.dist(&copy_nums_true),
                p.to_log_value(),
                format!(
                    "({:<3}{:<3}),({:<3}{:<3})",
                    n_missing, n_missing_null, n_error, n_error_null,
                ),
                dbg.summarize_cycle_with_dir(cycle),
            ),
            dbg,
        );
    });

    // (a)
    let neighbors: Vec<_> = neighbors
        .into_iter()
        .map(|(copy_nums, _, _, p)| (copy_nums, p))
        .collect();
    dbg_true.inspect_kmer_variance(&neighbors, &copy_nums_true);
}

fn main() {
    run_mcmc();
}
