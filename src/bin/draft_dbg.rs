use clap::Clap;
use dbgphmm::{
    dbg::{Dbg, HashDbg, SimpleDbg},
    e2e::{generate_dataset, Dataset, ReadType},
    genome,
    kmer::VecKmer,
    prelude::*,
};
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

#[derive(Clap)]
struct Opts {
    dataset_json: PathBuf,
    #[clap(short = 'k', default_value = "40")]
    k: usize,
}

fn main() {
    let opts: Opts = Opts::parse();
    let file = File::open(&opts.dataset_json).unwrap();
    let reader = BufReader::new(file);
    let dataset: Dataset = serde_json::from_reader(reader).unwrap();
    eprintln!("constructing..");
    let mut dbg_raw: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(opts.k, dataset.reads());
    dbg_raw.remove_nodes(2);
    dbg_raw.remove_deadend_nodes();
    println!("n_nodes={}", dbg_raw.n_nodes());
    println!("n_edges={}", dbg_raw.n_edges());
    println!("{:?}", dbg_raw.copy_num_stats());
    println!("{:?}", dbg_raw.degree_stats());

    // assert that true-kmers in the graph
    let (copy_nums_true, _) = dbg_raw
        .to_copy_nums_of_styled_seqs(dataset.genome())
        .unwrap();
    dbg_raw.set_node_copy_nums(&copy_nums_true);

    /*
    eprintln!("shrinking..");
    let mut dbg_shrinked = dbg_raw.shrink_single_copy_node();
    eprintln!("removing..");
    dbg_shrinked.remove_zero_copy_node();
    eprintln!("dumping..");
    println!("{}", dbg_shrinked.to_json());
    */
}
