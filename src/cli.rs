use crate::dbg::DBG;
use crate::hmm::base::PHMM;
use crate::hmm::params::PHMMParams;
use crate::hmm::sampler::PHMMSampler;
use crate::kmer::kmer::Kmer;
use crate::*;
use log::{info, warn};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;

pub fn generate(length: usize, seed: u64) {
    let v = random_seq::generate(length, seed);
    println!(">randseq");
    println!("{}", std::str::from_utf8(&v).unwrap());
}

pub fn sample(dbg_fa: String, length: u32, n_reads: u32, k: usize, seed: u64, param: PHMMParams) {
    let seqs = io::fasta::parse_seqs(&dbg_fa);
    let d = hmm::dbg::DbgPHMM::from_seqs(seqs, k);
    // println!("{}", d.as_dot());
    // println!("{}", d.dbg.as_dot());
    // let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    for i in 0..n_reads {
        // let seed_for_a_read: u64 = rng.gen();
        let seed_for_a_read = seed + i as u64;
        let seq = d.sample(&param, length, seed_for_a_read);
        let id = format!("{},{}", length, seed_for_a_read);
        io::fasta::dump_seq(&id, &seq);
    }
}

pub fn calc_prob(dbg_fa: String, reads_fa: String, k: usize, param: PHMMParams) {
    let seqs = io::fasta::parse_seqs(&dbg_fa);
    let d = hmm::dbg::DbgPHMM::from_seqs(seqs, k);
    let reads = io::fasta::parse_seqs(&reads_fa);

    let p = d.forward_prob(&param, &reads[0]);
    println!("forward prob : {}", p);
    // let p = d.backward(&param, &reads[0]);
    // println!("backward prob : {}", p[0].FM);
}

pub fn optimize(dbg_fa: String, reads_fa: String, k: usize, param: PHMMParams) {
    let seqs = io::fasta::parse_seqs(&dbg_fa);
    let d = hmm::dbg::DbgPHMM::from_seqs(seqs, k);
    info!("{}", d.dbg.as_degree_stats());

    let reads = io::fasta::parse_seqs(&reads_fa);
    let d = hmm::dbg::DbgPHMM::from_seqs(reads, k);
    info!("{}", d.dbg.as_degree_stats());

    let root = kmer::kmer::null_kmer(k - 1);
    info!("root={}", root);

    let s = cycles::DbgTree::new(&d.dbg, &root);
    info!("{}", s.as_stats());
}

pub fn sandbox() {
    let mut d = dbg::DbgHash::new();
    /*
    let seq = b"ATCGATTCGATTCGAT";
    let d = dbg::DbgHash::from_seq(seq, 5);
    */
    d.add_seq(b"ATCGATTCGATTCGAT", 8);
    // d.add_seq(b"ATCGATGCGATTCGAT", 8);
    println!("{}", d.as_dot());
    // eprintln!("{}", d.as_degree_stats());
    // eprintln!("{}", d.is_copy_number_consistent());
    let root = Kmer::from(b"NNNNNNNA");
    let s = cycles::DbgTree::new(&d, &root);
    for e in s.cycle_keys().iter() {
        println!("cycle {} = {:?}", e, s.cycle_components(e));
    }
}

pub fn sandbox2() {
    optimizer::base::test();
}

pub fn sandbox3() {
    // compressed_dbg::test();
    // optimizer::cdbg::test();
    hmm::cdbg::test();
}
