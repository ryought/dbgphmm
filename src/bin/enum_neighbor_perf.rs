use clap::{AppSettings, ArgEnum, Clap};
use dbgphmm::e2e::{generate_dataset, Dataset, ReadType};
use dbgphmm::genome;
use dbgphmm::graph::cycle::CycleWithDir;
use dbgphmm::kmer::common::kmers_to_string;
use dbgphmm::prelude::*;
use rayon::prelude::*;
use std::time::{Duration, Instant};

fn main() {
    let (genome, genome_size) = genome::tandem_repeat_haploid(20, 100, 0.01, 0, 0);
    let coverage = 10;
    let param = PHMMParams::uniform(0.001);
    let dataset = generate_dataset(
        genome.clone(),
        genome_size,
        0,
        param,
        coverage,
        5000,
        ReadType::FullLength,
        16,
        40,
    );
    let dbg_raw = dataset.dbg_raw.clone();
    // let mut dbg_true = dbg_raw.clone();
    let mut dbg_true = dbg_raw.clone().shrink_single_copy_node();
    dbg_true.remove_zero_copy_node();
    let (copy_nums_true, _) = dbg_true
        .to_copy_nums_of_styled_seqs(&genome)
        .unwrap_or_else(|err| panic!("true kmers {} is missing", kmers_to_string(&err)));
    dbg_true.set_node_copy_nums(&copy_nums_true);

    println!("# genome={}", genome[0]);
    println!("# k={}", dbg_true.k());
    println!("# n_kmers_with_null={}", dbg_true.n_kmers_with_null());
    println!("# n_dead_nodes={}", dbg_true.n_dead_nodes());
    println!("# n_nodes={}", dbg_true.n_nodes());
    println!("# n_edges={}", dbg_true.n_edges());

    let start = Instant::now();
    let neighbors = dbg_true.neighbor_copy_nums_fast();
    let duration = start.elapsed();
    println!("# n_neighbors={}", neighbors.len());
    println!("# time_neighbors={}", duration.as_millis());
    for neighbor in neighbors {
        println!("c={}", neighbor.sum());
    }

    // let start = Instant::now();
    // let neighbors = dbg_true.neighbor_copy_nums();
    // let duration = start.elapsed();
    // println!("# n_neighbors={}", neighbors.len());
    // println!("# time_neighbors={}", duration.as_millis());
}
