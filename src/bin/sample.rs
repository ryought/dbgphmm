use clap::Parser;
use dbgphmm::{
    common::ei,
    e2e::Dataset,
    hmmv2::params::PHMMParams,
    multi_dbg::posterior::test::test_inference_from_dbg,
    multi_dbg::MultiDbg,
    utils::{check_memory_usage, timer},
};

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(short = 'p')]
    p_infer: f64,
    #[clap(short = 's', default_value = "200")]
    sigma: usize,
    #[clap(long)]
    dbg: std::path::PathBuf,
    #[clap(long)]
    dataset_json: std::path::PathBuf,
    #[clap(long, use_value_delimiter = true, value_delimiter = ',')]
    edges: Vec<usize>,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# n_threads={}", rayon::current_num_threads());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);

    let dataset = Dataset::from_json_file(opts.dataset_json);
    let dbg = MultiDbg::from_dbg_file(&opts.dbg);
    // let neighbors = dbg.to_rescue_neighbors(2, 10);
    let neighbors: Vec<_> = opts
        .edges
        .iter()
        .flat_map(|&i| {
            let edge = ei(i);
            [
                dbg.to_rescue_neighbors_for_edge(edge, 2, true, false),
                dbg.to_rescue_neighbors_for_edge(edge, 2, true, true),
                dbg.to_rescue_neighbors_for_edge(edge, 2, false, false),
            ]
            .concat()
        })
        .collect();
    println!("neighbors={}", neighbors.len());
    let posterior = dbg.sample_posterior_for_inspect(
        neighbors,
        PHMMParams::uniform(opts.p_infer),
        dataset.reads(),
        dataset.genome_size(),
        opts.sigma,
    );
    let paths_true = dbg.paths_from_styled_seqs(dataset.genome()).unwrap();
    let copy_nums_true = dbg.copy_nums_from_full_path(paths_true);
    dbg.to_inspect_file(
        opts.dbg.with_extension("inspect"),
        &posterior,
        Some(&copy_nums_true),
    );

    println!("# finished_at={}", chrono::Local::now());
}