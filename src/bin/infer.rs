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
struct Opts {
    #[clap(short = 'k')]
    k_init: Option<usize>,
    #[clap(short = 'K')]
    k_max: usize,
    #[clap(short = 'p')]
    p_infer: f64,
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
        opts.output_prefix,
    );

    println!("# finished_at={}", chrono::Local::now());
}
