use crate::hmm::base::PHMM;
use crate::hmm::params::PHMMParams;
use crate::hmm::sampler::PHMMSampler;
use crate::*;
use log::{info, warn};

pub fn generate(length: usize, seed: u64) {
    let v = random_seq::generate(length, seed);
    println!(">randseq");
    println!("{}", std::str::from_utf8(&v).unwrap());
}

pub fn sample(dbg_fa: String, length: u32, n_reads: u32, k: usize, param: PHMMParams) {
    println!("sampling {} {} {} {}", dbg_fa, length, n_reads, k);
    println!("with param {}", param);

    let (kmers, copy_nums) = io::fasta::parse_kmers_and_copy_nums(&dbg_fa, k);
    let d = hmm::dbg::DbgPHMM::new(kmers, copy_nums).unwrap();
    let es = d.sample(&param, 10, 0);
}

pub fn calc_prob(dbg_fa: String, reads_fa: String, k: usize, param: PHMMParams) {
    let (kmers, copy_nums) = io::fasta::parse_kmers_and_copy_nums(&dbg_fa, k);
    let reads = io::fasta::parse_reads(&reads_fa);

    info!("from dbg_fa #kmer:{}", kmers.len());
    let d = hmm::dbg::DbgPHMM::new(kmers, copy_nums).unwrap();
    let p = d.forward_prob(&param, &reads[0]);
    println!("forward prob : {}", p);
    // let p = d.backward(&param, &reads[0]);
    // println!("backward prob : {}", p[0].FM);
}

pub fn sandbox() {
    println!("hoge");
}
