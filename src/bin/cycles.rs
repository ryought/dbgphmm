use clap::Parser;
use dbgphmm::{
    e2e::Dataset,
    genome,
    hmmv2::params::PHMMParams,
    multi_dbg::posterior::test::test_posterior_from_true,
    multi_dbg::{MultiDbg, NeighborConfig},
    utils::{check_memory_usage, timer},
};
use petgraph::algo::connected_components;

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(long)]
    dbg: std::path::PathBuf,
    #[clap(long, default_value = "10")]
    k_shortest: usize,
    #[clap(long)]
    prohibit_zero_copy: bool,
    #[clap(long)]
    exclude_passing_terminal: bool,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# n_threads={}", rayon::current_num_threads());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);
    // let dataset = Dataset::from_json_file(&opts.dataset);
    let dbg = MultiDbg::from_dbg_file(&opts.dbg);

    for (edge, _, _, _) in dbg.edges_compact() {
        if dbg.copy_num_of_edge_in_compact(edge) == 0 {
            let neighbor_copy_nums: Vec<_> = dbg
                .to_rescue_neighbors_for_edge(edge, opts.k_shortest, opts.prohibit_zero_copy)
                .into_iter()
                .filter(|(_, info)| {
                    !opts.exclude_passing_terminal || !dbg.is_passing_terminal(&info)
                })
                .collect();
            println!("e{} n_cycles={}", edge.index(), neighbor_copy_nums.len());
            for (c, i) in neighbor_copy_nums {
                println!("\t{:?}", i);
            }
        }
    }
}
