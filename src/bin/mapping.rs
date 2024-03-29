//!
//! Show mapping of forward/backward on PHMM corresponding to the MultiDbg
//!
//! generate bmap file
//!
use clap::Parser;
use dbgphmm::{e2e::Dataset, hmmv2::params::PHMMParams, multi_dbg::MultiDbg, utils::timer};
use petgraph::graph::NodeIndex;

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(long)]
    dbg: std::path::PathBuf,
    #[clap(long)]
    dataset: std::path::PathBuf,
    #[clap(long)]
    map_input: Option<std::path::PathBuf>,
    #[clap(long)]
    map_output: Option<std::path::PathBuf>,
    #[clap(short = 'e')]
    p_error: f64,
    #[clap(short = 'r')]
    max_ratio: Option<f64>,
    #[clap(short = 'a')]
    n_active_nodes: usize,
    #[clap(long, arg_enum, default_value = "uniform")]
    phmm: PHMM,
}

#[derive(Clone, Debug, clap::ArgEnum)]
enum PHMM {
    Normal,
    Uniform,
    NonZero,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# n_threads={}", rayon::current_num_threads());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);
    let dataset = Dataset::from_json_file(&opts.dataset);
    let dbg = MultiDbg::from_dbg_file(&opts.dbg);

    let mut node_name_vec = vec![String::new(); dbg.n_edges_full()];
    for ec in dbg.graph_compact().edge_indices() {
        for (i, ef) in dbg.edges_in_full(ec).iter().enumerate() {
            let base = dbg.base(*ef);
            node_name_vec[ef.index()] = format!("{}-{}{}", ec.index(), i, base as char);
        }
    }
    let format = |node: NodeIndex| node_name_vec[node.index()].clone();

    let mut param = PHMMParams::uniform(opts.p_error);
    param.n_active_nodes = opts.n_active_nodes;
    param.active_node_max_ratio = opts.max_ratio.unwrap_or(30.0);

    let phmm = match opts.phmm {
        PHMM::Normal => dbg.to_phmm(param),
        PHMM::Uniform => dbg.to_uniform_phmm(param),
        PHMM::NonZero => dbg.to_non_zero_phmm(param),
    };

    // n_zero_edges
    let n_zero_edges = dbg
        .graph_compact()
        .edge_indices()
        .filter(|&e| dbg.copy_num_of_edge_in_compact(e) == 0)
        .count();
    println!("# n_zero_edges={}", n_zero_edges);

    if let Some(map) = &opts.map_output {
        let mappings = dbg.generate_mappings(param, dataset.reads(), None);
        let freqs = dbg.mappings_to_freqs(&mappings);
        let copy_num =
            dbg.min_squared_error_copy_nums_from_freqs(&freqs, dataset.coverage(), Some(2));
        eprintln!("copy_num={}", copy_num);
        dbg.to_map_file(map, dataset.reads(), &mappings).unwrap();
        eprintln!("map {} written", map.display());
        return;
    }

    let mappings = if let Some(map_filename) = &opts.map_input {
        Some(dbg.from_map_file(map_filename, dataset.reads()))
    } else {
        None
    };

    for (i, read) in dataset.reads().into_iter().enumerate() {
        let (output, t) = timer(|| {
            if mappings.is_some() {
                println!("# using mapping");
                phmm.run_with_mapping(read, &mappings.as_ref().unwrap()[i])
            } else if opts.max_ratio.is_some() {
                println!("# using max_ratio");
                phmm.run_sparse_adaptive(read, true)
            } else {
                println!("# using n_active_nodes");
                phmm.run_sparse_adaptive(read, false)
            }
        });

        // summary
        println!(
            "# i={} L={} t={} pF={} pB={}",
            i,
            read.len(),
            t,
            output.to_full_prob_forward().to_log_value(),
            output.to_full_prob_backward().to_log_value(),
        );

        // for each base
        for j in 0..output.n_emissions() {
            let base = read.seq()[j] as char;
            println!("# {}\t{}\tR\t{}", i, j, base);
            let f = output.forward.table_merged(j + 1);
            println!(
                "{}\t{}\tF\t{}\t{}",
                i,
                j,
                f.e.to_log_value(),
                f.to_summary_string_n(40, format),
            );
            let b = output.backward.table_merged(j + 1);
            println!(
                "{}\t{}\tB\t{}\t{}",
                i,
                j,
                b.mb.to_log_value(),
                b.to_summary_string_n(40, format),
            );
            println!(
                "{}\t{}\tS\t{}\t{}",
                i,
                j,
                f.n_active_nodes(),
                output.to_emit_probs(j + 1).to_summary_string_n(40, format),
            );
        }
    }
    println!("# finished_at={}", chrono::Local::now());
}
