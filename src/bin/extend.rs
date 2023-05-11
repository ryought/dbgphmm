use clap::Parser;
use dbgphmm::{
    e2e::Dataset,
    genome,
    hmmv2::params::PHMMParams,
    multi_dbg::posterior::test::test_inference_from_dbg,
    multi_dbg::MultiDbg,
    utils::{check_memory_usage, timer},
};

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(short = 'k')]
    k_init: usize,
    #[clap(short = 'K')]
    k_max: usize,
    #[clap(long)]
    dataset_json: std::path::PathBuf,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# n_threads={}", rayon::current_num_threads());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);

    let dataset = Dataset::from_json_file(opts.dataset_json);
    let dbg = MultiDbg::create_draft_from_dataset(opts.k_init, &dataset);
    let mappings = dbg.generate_mappings(dataset.params(), dataset.reads(), None);
    dbg.purge_and_extend(&[], opts.k_max, false, None, &mappings);

    println!("# finished_at={}", chrono::Local::now());
}
