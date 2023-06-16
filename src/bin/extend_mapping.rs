use clap::Parser;
use dbgphmm::{
    e2e::Dataset, genome, hmmv2::params::PHMMParams,
    multi_dbg::posterior::test::test_mapping_extension, multi_dbg::MultiDbg,
};

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(short = 'K')]
    k_max: usize,
    #[clap(short = 'p')]
    p_infer: f64,
    #[clap(long)]
    dbg: std::path::PathBuf,
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
    let dbg = MultiDbg::from_dbg_file(opts.dbg);
    test_mapping_extension(
        &dataset,
        dbg,
        opts.k_max,
        PHMMParams::uniform(opts.p_infer),
        opts.output_prefix,
    );

    println!("# finished_at={}", chrono::Local::now());
}
