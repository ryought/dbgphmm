use clap::{AppSettings, Clap};
use dbgphmm::hmm::base::PHMM;
use dbgphmm::prob::Prob;
use dbgphmm::*;
use std::io::prelude::*;

fn test() {
    // generics
    let l1 = vec![34, 50, 25];
    println!("l1: {}", hoge::fuga::largest(&l1));

    // string
    hoge::mystring::test1();

    // DNA
    /*
    let x = hoge::seq::MyDNA::A;
    let c = hoge::seq::show(x);
    let c2 = hoge::seq::show(hoge::seq::complement(x));
    println!("c = {}, c2 = {}", c, c2);
    */

    test_struct::hoge::test1();
}

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
    let reads = io::fasta::parse_reads(&opts.reads_fa);
    let d = hmm::dbg::DbgPHMM::new(kmers, copy_nums).unwrap();
    let param = hmm::params::PHMMParams::new(
        Prob::from_prob(0.01),
        Prob::from_prob(0.01),
        Prob::from_prob(0.01),
        10,
    );
    // let p = d.forward_prob(&param, &reads[0]);
    // println!("forward prob : {}", p);
    // let es = d.sample(&param, 10);
    // println!("emmissions: {:?}", es);
    // hmm::base::test_random();
    hmm::testing::test_static();
}

fn run2(config: kmer::counter::Config) {
    let mut f = std::fs::File::open(config.filename).expect("file not found");
    let mut contents = String::new();
    f.read_to_string(&mut contents).expect("cannot read file");
    println!("text: {}", contents);
}

fn length(seq: &[u8]) -> usize {
    // seq.len()
    let mut n = 0;
    for s in seq.iter() {
        n += 1;
    }
    n
}

use std::collections::HashMap;
fn kmer_count(seq: &[u8], k: usize) -> HashMap<&[u8], usize> {
    let mut count = HashMap::new();
    for i in 0..10 {
        count.insert(&seq[i..i + k], 1);
    }
    count
}

#[derive(Debug)]
enum DNA {
    A,
    C,
    G,
    T,
    N,
}
#[derive(Debug)]
struct Seq {
    array: Vec<DNA>,
}
/*
impl Seq {
    fn from_slice() -> Seq {
    }
}
*/
fn test5() {
    let s1 = Seq {
        array: vec![DNA::A, DNA::C, DNA::T],
    };
    println!("seq {:?}", s1)
}
