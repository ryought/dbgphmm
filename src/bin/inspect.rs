use clap::Parser;
use dbgphmm::{
    e2e::Dataset,
    genome,
    hmmv2::params::PHMMParams,
    multi_dbg::posterior::test::test_inference_from_dbg,
    multi_dbg::{MultiDbg, NeighborConfig},
    utils::{check_memory_usage, timer},
};

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    dbg: std::path::PathBuf,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);

    let dbg = MultiDbg::from_dbg_file(opts.dbg);
    for (copy_nums, info) in dbg.to_neighbor_copy_nums_and_infos(NeighborConfig {
        max_cycle_size: 10,
        max_flip: 2,
        use_long_cycles: true,
        ignore_cycles_passing_terminal: false,
        use_reducers: true,
    }) {
        println!("{} {:?}", copy_nums, info);
    }
}
