use clap::{AppSettings, Clap};
use dbgphmm::hmm::base::PHMM;
use dbgphmm::prob::Prob;
use dbgphmm::*;
use log::{info, warn};
use std::io::prelude::*;

/// de bruijn graph + profile HMM optimization package
#[derive(Clap)]
#[clap(version = "0.1", author = "ryought <ryonakabayashi@gmail.com>")]
#[clap(setting = AppSettings::ColoredHelp)]
struct Opts {
    /// Fasta input
    #[clap(default_value = "test.fa")]
    dbg_fa: String,
    /// Read fasta input
    #[clap(default_value = "read.fa")]
    reads_fa: String,
    /// initial kmer size
    #[clap(short, default_value = "8")]
    k: usize,
    /// Print debug info
    #[clap(short, long)]
    debug: bool,
}

fn main() {
    // enable logger
    env_logger::init();

    // parse options
    let opts: Opts = Opts::parse();
    let (kmers, copy_nums) = io::fasta::parse_kmers_and_copy_nums(&opts.dbg_fa, opts.k);
    info!("from dbg_fa #kmer:{}", kmers.len());
    let d = hmm::dbg::DbgPHMM::new(kmers, copy_nums).unwrap();
    let reads = io::fasta::parse_reads(&opts.reads_fa);
    let param = hmm::params::PHMMParams::new(
        Prob::from_prob(0.01),
        Prob::from_prob(0.01),
        Prob::from_prob(0.01),
        3,
    );
    let p = d.forward_prob(&param, &reads[0]);
    println!("forward prob : {}", p);
    // let es = d.sample(&param, 10);
    // println!("emmissions: {:?}", es);
    // hmm::base::test_random();
    // hmm::testing::test_static();
    // let v = random_seq::generate(100, 0);
    // println!("{}", std::str::from_utf8(&v).unwrap());
}
