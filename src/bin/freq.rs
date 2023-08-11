use clap::Parser;
use dbgphmm::{
    e2e::Dataset,
    genome,
    hmmv2::params::PHMMParams,
    multi_dbg::posterior::test::test_inference_from_dbg,
    multi_dbg::MultiDbg,
    utils::{check_memory_usage, timer},
};
use petgraph::graph::NodeIndex;

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(short = 'e')]
    p_error: f64,
    #[clap(long)]
    dbg: std::path::PathBuf,
    #[clap(long)]
    dataset_json: std::path::PathBuf,
    // #[clap(long)]
    // output_prefix: std::path::PathBuf,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# n_threads={}", rayon::current_num_threads());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);

    let dataset = Dataset::from_json_file(opts.dataset_json);
    let mut dbg = MultiDbg::from_dbg_file(&opts.dbg);
    let param = PHMMParams::uniform(opts.p_error);

    eprintln!("mapping...");
    let (mappings, t) = timer(|| dbg.generate_mappings(param, dataset.reads(), None));
    println!("# map t={}", t);
    let freqs = dbg.mappings_to_freqs(&mappings);
    eprintln!("approximate...");
    let (copy_nums, t) =
        timer(|| dbg.min_squared_error_copy_nums_from_freqs(&freqs, dataset.coverage(), Some(2)));
    dbg.set_copy_nums(&copy_nums);
    println!("# approx t={}", t);
    println!("# {}", copy_nums);

    eprintln!("calculate true...");
    let paths_true = dbg.paths_from_styled_seqs(dataset.genome()).unwrap();
    let copy_nums_true = dbg.copy_nums_from_full_path(paths_true);
    println!("# {}", copy_nums_true);

    let mut n_over_estimate = 0;
    let mut n_under_estimate = 0;
    let mut sum_over_estimate = 0;
    let mut sum_under_estimate = 0;

    for (i, e) in dbg
        .graph_compact()
        .edge_indices()
        .filter(|&e| copy_nums[e] != copy_nums_true[e])
        .enumerate()
    {
        let freq_sum: f64 = dbg
            .edges_in_full(e)
            .iter()
            .map(|e| freqs[NodeIndex::new(e.index())])
            .sum();
        let freq_ave = freq_sum / (dbg.edges_in_full(e).len() as f64);

        println!(
            "{} e{} {} {} {}",
            i,
            e.index(),
            copy_nums[e],
            copy_nums_true[e],
            freq_ave / dataset.coverage()
        );

        if copy_nums[e] > copy_nums_true[e] {
            n_over_estimate += 1;
            sum_over_estimate += copy_nums_true[e];
        } else {
            n_under_estimate += 1;
            sum_under_estimate += copy_nums_true[e];
        }
    }

    println!("over n={} sum={}", n_over_estimate, sum_over_estimate);
    println!("under n={} sum={}", n_under_estimate, sum_under_estimate);

    println!("# finished_at={}", chrono::Local::now());
}
