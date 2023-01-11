use clap::Clap;
use dbgphmm::hmmv2::params::PHMMParams;
use dbgphmm::inspect::{generate_500bp_case, generate_small_case};
use std::fs::File;
use std::io::prelude::*;
use std::path::PathBuf;

#[derive(Clap, Debug)]
struct Opts {
    #[clap(short)]
    a: usize,
    #[clap(short)]
    b: usize,
    #[clap(short)]
    c: usize,
    #[clap(short)]
    d: usize,
    #[clap(short)]
    p: f64,
    #[clap(long)]
    pi: f64,
    #[clap(long)]
    dbgviz_output: Option<PathBuf>,
    #[clap(long)]
    show_mapping: bool,
}

fn main() {
    let opts: Opts = Opts::parse();
    let (dataset, dbg_true, dbg_opt) = generate_small_case(opts.a, opts.b, opts.c, opts.d, opts.p);
    dataset.show_genome();
    // dataset.show_reads_with_genome();
    println!("{}", dbg_true);
    println!("{}", dbg_opt);
    let copy_nums_true = dbg_true.to_node_copy_nums();
    let copy_nums_opt = dbg_opt.to_node_copy_nums();

    let param = PHMMParams::uniform(opts.pi);
    let p_true = dbg_true.to_full_prob(param, dataset.genome());
    let p_opt = dbg_opt.to_full_prob(param, dataset.genome());
    println!("p_true(g)={}", p_true);
    println!("p_opt(g)={}", p_opt);

    for i in 0..(dataset.genome().len()) {
        let p_true = dbg_true.to_full_prob(param, &[&dataset.genome()[i]]);
        println!("p_true(g[{}])={}", i, p_true);
        let p_opt = dbg_opt.to_full_prob(param, &[&dataset.genome()[i]]);
        println!("p_opt(g[{}])={}", i, p_opt);
    }

    if opts.show_mapping {
        for i in 0..(dataset.genome().len()) {
            let p_true = dbg_true.to_full_prob(param, &[&dataset.genome()[i]]);
            println!("p_true(g[{}])={}", i, p_true);
            let p_opt = dbg_opt.to_full_prob(param, &[&dataset.genome()[i]]);
            println!("p_opt(g[{}])={}", i, p_opt);

            dbg_true.compare_mappings_v2(
                param,
                dataset.genome()[i].seq(),
                &copy_nums_true,
                &copy_nums_opt,
            );
        }
    }

    if let Some(path) = &opts.dbgviz_output {
        let json = dbg_true.to_cytoscape_with_info(
            |node| {
                Some(format!(
                    "{}x {}x",
                    copy_nums_true[node], copy_nums_opt[node]
                ))
            },
            None,
        );
        let mut file = File::create(path).unwrap();
        writeln!(file, "{}", json).unwrap();
    }
}
