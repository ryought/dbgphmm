use clap::Parser;
use dbgphmm::{
    e2e::Dataset,
    genome,
    hmmv2::params::PHMMParams,
    multi_dbg::{MultiDbg, NeighborConfig},
    utils::{check_memory_usage, timer},
};
use petgraph::algo::connected_components;

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(long)]
    dbg: std::path::PathBuf,
    #[clap(long)]
    dataset: std::path::PathBuf,
    #[clap(long)]
    map: std::path::PathBuf,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# n_threads={}", rayon::current_num_threads());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);

    let dataset = Dataset::from_json_file(opts.dataset);
    let dbg = MultiDbg::from_dbg_file(opts.dbg);
    println!("k={} |E|={}", dbg.k(), dbg.n_edges_full());
    let mappings = dbg.from_map_file(opts.map, dataset.reads());

    let mut params = dataset.params();
    params.n_warmup = dbg.k() + 5;
    println!("params={}", params);
    let phmm = dbg.to_phmm(params);
    // read index
    let r = 0;

    let o_with_map = phmm.run_with_mapping(&dataset.reads()[r], &mappings[r]);
    println!(
        "with_map F={} B={}",
        o_with_map.to_full_prob_forward(),
        o_with_map.to_full_prob_backward()
    );

    let o_without_map = phmm.run_sparse(&dataset.reads()[r]);
    println!(
        "without_map F={} B={}",
        o_without_map.to_full_prob_forward(),
        o_without_map.to_full_prob_backward()
    );

    for i in 0..o_with_map.n_emissions() {
        println!("i={}", i);
        println!(
            "with\t{}",
            o_with_map
                .forward
                .table(i)
                .to_summary_string_n(100, |_| String::new())
        );
        println!(
            "without\t{}",
            o_without_map
                .forward
                .table(i)
                .to_summary_string_n(100, |_| String::new())
        );
        println!("pF_with={}", o_with_map.forward.table(i).e);
        println!("pF_without={}", o_without_map.forward.table(i).e);
        println!("map\t{}", mappings[r].to_nodes_string(i));

        println!(
            "with\t{}",
            o_with_map
                .backward
                .table(i)
                .to_summary_string_n(100, |_| String::new())
        );
        println!(
            "without\t{}",
            o_without_map
                .backward
                .table(i)
                .to_summary_string_n(100, |_| String::new())
        );
        println!("pB_with={}", o_with_map.backward.table(i).mb);
        println!("pB_without={}", o_without_map.backward.table(i).mb);
    }
}
