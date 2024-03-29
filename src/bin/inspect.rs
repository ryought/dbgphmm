use clap::Parser;
use dbgphmm::{e2e::Dataset, multi_dbg::MultiDbg};

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    dbg: std::path::PathBuf,
    dataset_json: std::path::PathBuf,
    #[clap(short = 's', default_value = "200")]
    sigma: usize,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# n_threads={}", rayon::current_num_threads());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);

    let dataset = Dataset::from_json_file(opts.dataset_json);
    let mut dbg = MultiDbg::from_dbg_file(opts.dbg);
    println!("k={} |E|={}", dbg.k(), dbg.n_edges_full());
    let mut params = dataset.params();
    params.n_warmup = dbg.k() + 5;
    println!("params={}", params);

    // original
    let score = dbg.to_score(
        params,
        dataset.reads(),
        None,
        dataset.genome_size(),
        opts.sigma,
    );
    let likelihood = dbg.to_phmm(params).to_full_prob_parallel(dataset.reads());
    println!(
        "orig\t{}\t{}\t{}\t{}",
        score.p(),
        likelihood,
        score,
        dbg.get_copy_nums()
    );
    // println!("orig\tmapping...");
    // let mappings = dbg.generate_mappings(params, dataset.reads(), None);
    // dbg.to_map_file("orig.map", dataset.reads(), &mappings);

    // println!("orig\tncc={}", connected_components(dbg.graph_compact()));

    // true
    let paths_true = dbg
        .paths_from_styled_seqs(dataset.genome())
        .expect("k-mer in genome is missing");
    let copy_nums_true = dbg.copy_nums_from_full_path(paths_true);
    dbg.set_copy_nums(&copy_nums_true);
    let score = dbg.to_score(
        params,
        dataset.reads(),
        None,
        dataset.genome_size(),
        opts.sigma,
    );
    let likelihood = dbg.to_phmm(params).to_full_prob_parallel(dataset.reads());
    println!(
        "true\t{}\t{}\t{}\t{}",
        score.p(),
        likelihood,
        score,
        copy_nums_true
    );
    println!("true\tmapping...");
    let mappings = dbg.generate_mappings(params, dataset.reads(), None);
    dbg.to_map_file("true.map", dataset.reads(), &mappings)
        .unwrap();
    // println!("true\tncc={}", connected_components(dbg.graph_compact()));

    // for (copy_nums, info) in dbg.to_neighbor_copy_nums_and_infos(NeighborConfig {
    //     max_cycle_size: 10,
    //     max_flip: 2,
    //     use_long_cycles: true,
    //     ignore_cycles_passing_terminal: false,
    //     use_reducers: true,
    // }) {
    //     println!("{} {:?}", copy_nums, info);
    // }
}
