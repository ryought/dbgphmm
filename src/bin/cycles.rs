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
                .to_rescue_neighbors_for_edge(edge, 5, false)
                .into_iter()
                .filter(|(_, info)| !dbg.is_passing_terminal(&info))
                .filter(|(_, info)| dbg.has_zero_to_one_change(&info))
                .collect();
            println!("e{} n_cycles={}", edge.index(), neighbor_copy_nums.len());
            for (c, i) in neighbor_copy_nums {
                println!("{:?}", i);
            }
        }
    }
}
