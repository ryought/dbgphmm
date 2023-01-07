use clap::Clap;
use dbgphmm::inspect::{generate_500bp_case, generate_small_case};
use std::fs::File;
use std::io::prelude::*;

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
}

fn main() {
    let opts: Opts = Opts::parse();
    let (dataset, dbg_true, dbg_opt) = generate_small_case(opts.a, opts.b, opts.c, opts.d);
    // dataset.show_reads_with_genome();
    println!("{}", dbg_true);
    println!("{}", dbg_opt);
    let copy_nums_true = dbg_true.to_node_copy_nums();
    let copy_nums_opt = dbg_opt.to_node_copy_nums();
    let p_true = dbg_true.to_full_prob(dataset.params(), dataset.reads());
    println!("p_true={}", p_true);
    let p_opt = dbg_opt.to_full_prob(dataset.params(), dataset.reads());
    println!("p_opt={}", p_opt);

    // let json = dbg_true.to_cytoscape_with_info(
    //     |node| {
    //         Some(format!(
    //             "{}x {}x",
    //             copy_nums_true[node], copy_nums_opt[node]
    //         ))
    //     },
    //     None,
    // );
    // let mut file = File::create("inspect500bp_true_k12.json").unwrap();
    // writeln!(file, "{}", json).unwrap();
}
