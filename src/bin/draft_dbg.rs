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
    #[clap(long)]
    use_true_copy_nums: bool,
}

fn main() {
    let opts: Opts = Opts::parse();
    let file = File::open(&opts.dataset_json).unwrap();
    let reader = BufReader::new(file);
    let dataset: Dataset = serde_json::from_reader(reader).unwrap();
    eprintln!("constructing..");
    let mut dbg_draft: SimpleDbg<VecKmer> =
        SimpleDbg::create_draft_from_seqs(opts.k, dataset.reads(), dataset.coverage());

    // assert that true-kmers in the graph
    if opts.use_true_copy_nums {
        eprintln!("assigning true copynums..");
        dbg_draft.set_copy_nums_by_styled_seq(dataset.genome());
    }

    // output
    eprintln!("dumping..");
    println!("{}", dbg_draft.to_json());
}
