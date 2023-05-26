//!
//! Show mapping of forward/backward on PHMM corresponding to the MultiDbg
//!
use clap::Parser;
use dbgphmm::{
    e2e::Dataset,
    genome,
    hmmv2::params::PHMMParams,
    multi_dbg::posterior::test::test_posterior_from_true,
    multi_dbg::{MultiDbg, NeighborConfig},
};

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(long)]
    dbg: std::path::PathBuf,
    #[clap(long)]
    dataset: std::path::PathBuf,
    #[clap(short = 'e')]
    p_error: f64,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# n_threads={}", rayon::current_num_threads());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);
    let dataset = Dataset::from_json_file(&opts.dataset);
    let dbg = MultiDbg::from_dbg_file(&opts.dbg);

    let param = PHMMParams::uniform(opts.p_error);
    let mappings = dbg.generate_mappings(param, dataset.reads(), None);
    // let freqs = dbg.mappings_to_freqs(&mappings);
    // let copy_num = dbg.min_squared_error_copy_nums_from_freqs(&freqs, dataset.coverage(), Some(2));
    // println!("copy_num={}", copy_num);
}
