use clap::Parser;
use dbgphmm::{
    e2e::Dataset, hmmv2::params::PHMMParams, multi_dbg::posterior::test::test_posterior_from_true,
    multi_dbg::MultiDbg,
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
    println!("k={} |E|={}", dbg.k(), dbg.n_edges_full());
    let (mut dbg, post, copy_nums_true) = test_posterior_from_true(
        &dataset,
        dbg,
        500,
        opts.dbg.with_extension("from_true.inspect"),
        PHMMParams::uniform(opts.p_error),
    );
    dbg.set_copy_nums(&post.max_copy_nums());
    dbg.to_gfa_post_file(
        opts.dbg.with_extension("from_true.gfa"),
        &post,
        Some(&copy_nums_true),
    )
    .unwrap();
}
