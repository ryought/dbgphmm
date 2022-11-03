use dbgphmm::e2e::{generate_dataset, Dataset, ReadType};
use dbgphmm::genome;
use dbgphmm::graph::cycle::CycleWithDir;
use dbgphmm::prelude::*;
use rayon::prelude::*;

fn run_mcmc() {
    let (genome, genome_size) = genome::simple(100, 5);
    // let (genome, genome_size) = genome::tandem_repeat_haploid(20, 5, 0.01, 0, 0);
    // let (genome, genome_size) = genome::tandem_repeat_haploid(50, 4, 0.05, 0, 0);
    let coverage = 20;
    let param = PHMMParams::uniform(0.003);
    let dataset = generate_dataset(
        genome.clone(),
        genome_size,
        0,
        param,
        coverage,
        2000,
        ReadType::FullLength,
        32,
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

    let neighbors: Vec<_> = neighbors
        .into_par_iter()
        .filter(|(copy_nums, _)| copy_nums.sum() > 0)
        .map(|(copy_nums, cycle)| {
            let mut dbg = dbg_true.clone();
            dbg.set_node_copy_nums(&copy_nums);
            let p = dbg.to_full_prob(param, &dataset.reads);
            let ((n_missing, n_missing_null), (n_error, n_error_null)) = dbg.inspect_kmers(&genome);
            println!(
                "N\t{}\t{}\t{}\t({},{}),({},{})\t{}\t{}\t{}",
                dbg.genome_size(),
                copy_nums.dist(&copy_nums_true),
                p.to_log_value(),
                n_missing,
                n_missing_null,
                n_error,
                n_error_null,
                dbg.summarize_cycle_with_dir(&cycle),
                dbg,
                copy_nums,
            );
            (copy_nums, p)
        })
        .collect();

    // (a)
    dbg_true.inspect_kmer_variance(&neighbors);
}

fn main() {
    run_mcmc();
}
