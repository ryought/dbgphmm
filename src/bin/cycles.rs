use clap::Parser;
use dbgphmm::{
    e2e::Dataset,
    genome,
    hmmv2::params::PHMMParams,
    multi_dbg::posterior::test::test_posterior_from_true,
    multi_dbg::{MultiDbg, NeighborConfig},
    utils::{check_memory_usage, timer},
};
use petgraph::algo::connected_components;

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(long)]
    dbg: std::path::PathBuf,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# n_threads={}", rayon::current_num_threads());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);
    // let dataset = Dataset::from_json_file(&opts.dataset);
    let dbg = MultiDbg::from_dbg_file(&opts.dbg);

    let neighbor_copy_nums: Vec<_> = dbg
        .to_neighbor_copy_nums_and_infos(NeighborConfig {
            max_cycle_size: 50,
            max_flip: 4,
            use_long_cycles: true,
            ignore_cycles_passing_terminal: true,
            use_reducers: false,
        })
        .into_iter()
        .filter(|(_, info)| !dbg.is_passing_terminal(&info))
        .filter(|(_, info)| dbg.has_zero_to_one_change(&info))
        .collect();
    println!("n_cycles={}", neighbor_copy_nums.len());
}
