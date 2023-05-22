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
    k_init: Option<usize>,
    #[clap(short = 'K')]
    k_max: usize,
    #[clap(short = 'e')]
    p_error: f64,
    #[clap(short = 'p')]
    p_infer: f64,
    #[clap(short = 's', default_value = "200")]
    sigma: usize,
    #[clap(short = 'I', default_value = "10")]
    max_iter: usize,
    #[clap(short = 'c', default_value = "10")]
    max_cycle_size: usize,
    #[clap(long = "kr", default_value = "128")]
    k_max_rescue_only: usize,
    #[clap(long = "km", default_value = "100")]
    k_max_rerun_mapping: usize,
    #[clap(long)]
    dbg: Option<std::path::PathBuf>,
    #[clap(long)]
    dataset_json: std::path::PathBuf,
    #[clap(long)]
    output_prefix: std::path::PathBuf,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# n_threads={}", rayon::current_num_threads());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);

    let dataset = Dataset::from_json_file(opts.dataset_json);
    let dbg = if let Some(dbg_filename) = opts.dbg {
        MultiDbg::from_dbg_file(dbg_filename)
    } else {
        MultiDbg::create_draft_from_dataset(opts.k_init.unwrap(), &dataset)
    };
    test_inference_from_dbg(
        &dataset,
        dbg,
        opts.k_max,
        PHMMParams::uniform(opts.p_infer),
        PHMMParams::uniform(opts.p_error),
        opts.sigma,
        opts.max_iter,
        opts.max_cycle_size,
        opts.output_prefix,
        opts.k_max_rescue_only,
        opts.k_max_rerun_mapping,
    );

    println!("# finished_at={}", chrono::Local::now());
}
