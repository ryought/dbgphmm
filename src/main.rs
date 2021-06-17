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

/*
fn test() {
    println!("run!");
    let mut x: f64 = 0.0;
    for i in 0..100_000 {
        if i % 2 == 0 {
            x += (i as f64).ln()
        } else {
            x -= (i as f64).ln()
        }
    }
    println!("finish!");
}
*/

fn head_ref(vec: &[u8], i: usize) -> &[u8] {
    if i == 0 {
        &[]
    } else {
        &vec[i..i + 1]
    }
}

fn test_head_ref() {
    let v: Vec<u8> = vec![10, 20, 30, 40, 50, 60];
    for i in 0..6 {
        let r = head_ref(&v, i);
        println!("head ref {}: {:?}", i, r);
    }
}

fn main() {
    // enable logger
    env_logger::init();

    // hmm::dbg::test();

    // ref_test::test();

    // parse options
    let opts: Opts = Opts::parse();
    let (kmers, copy_nums) = io::fasta::parse_kmers_and_copy_nums(&opts.dbg_fa, opts.k);
    info!("from dbg_fa #kmer:{}", kmers.len());
    let d = hmm::dbg::DbgPHMM::new(kmers, copy_nums).unwrap();
    let reads = io::fasta::parse_reads(&opts.reads_fa);
    let param = hmm::params::PHMMParams::default();
    let p = d.forward_prob(&param, &reads[0]);
    println!("forward prob : {}", p);

    // let es = d.sample(&param, 10);
    // println!("emmissions: {:?}", es);
    // hmm::base::test_random();
    // hmm::testing::test_static();
    // let v = random_seq::generate(100, 0);
    // println!("{}", std::str::from_utf8(&v).unwrap());
}
