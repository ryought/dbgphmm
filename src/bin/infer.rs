use clap::Parser;
use dbgphmm::{
    e2e::Dataset, hmmv2::params::PHMMParams,
    multi_dbg::posterior::test::test_inference_from_dbg_with_dataset, multi_dbg::MultiDbg,
    prob::Prob,
};

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(short = 'k')]
    k_init: Option<usize>,
    /// Minimum occurrence of k-mers in read
    #[clap(short = 'm', default_value_t = 2)]
    min_count: usize,
    /// Minimum occurrence of deadend k-mers in read
    #[clap(short = 'M')]
    min_deadend_count: Option<usize>,
    #[clap(short = 'K')]
    k_max: usize,
    #[clap(short = 'e')]
    p_error: f64,
    #[clap(short = 'p')]
    p_infer: f64,
    #[clap(long, default_value = "0.8")]
    p0: f64,
    #[clap(short = 's')]
    sigma: usize,
    #[clap(short = 'I', default_value = "10")]
    max_iter: usize,
    #[clap(short = 'c', default_value = "10")]
    max_cycle_size: usize,
    #[clap(long)]
    dbg: Option<std::path::PathBuf>,
    #[clap(long)]
    map: Option<std::path::PathBuf>,
    #[clap(long)]
    dataset_json: std::path::PathBuf,
    #[clap(long)]
    output_prefix: std::path::PathBuf,
    /// Number of threads used in parallel posterior sampling.
    /// Use all threads if not specified.
    #[clap(short = 't')]
    n_threads: Option<usize>,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);

    if let Some(n_threads) = opts.n_threads {
        println!("# n_threads_specified={}", n_threads);
        rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build_global()
            .unwrap();
    }
    println!("# n_threads={}", rayon::current_num_threads());

    let dataset = Dataset::from_json_file(opts.dataset_json);
    let dbg = if let Some(dbg_filename) = opts.dbg {
        MultiDbg::from_dbg_file(dbg_filename)
    } else {
        let min_deadend_count = opts
            .min_deadend_count
            .unwrap_or((dataset.coverage() / 4.0) as usize);
        MultiDbg::create_draft_from_dataset_with(
            opts.k_init.unwrap(),
            &dataset,
            opts.min_count,
            min_deadend_count,
        )
    };
    let mappings = if let Some(map_filename) = opts.map {
        Some(dbg.from_map_file(map_filename, dataset.reads()))
    } else {
        None
    };
    test_inference_from_dbg_with_dataset(
        &dataset,
        dbg,
        opts.k_max,
        PHMMParams::uniform(opts.p_infer),
        PHMMParams::uniform(opts.p_error),
        opts.sigma,
        opts.max_iter,
        opts.max_cycle_size,
        Prob::from_prob(opts.p0),
        opts.output_prefix,
        mappings,
    );

    println!("# finished_at={}", chrono::Local::now());
}
