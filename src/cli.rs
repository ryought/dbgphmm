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

pub fn stat(dbg_fa: String, k: usize) {
    let seqs = io::fasta::parse_seqs(&dbg_fa);
    let d = dbg::DbgHash::from_seqs(seqs, k);
    info!("deg {}", d.as_degree_stats());
    let cdbg = compressed_dbg::CompressedDBG::from(&d, k);
    let copy_nums_true: Vec<u32> = cdbg
        .iter_nodes()
        .map(|v| {
            let kmer = cdbg.kmer(&v);
            d.find(kmer)
        })
        .collect();

    let total = cdbg.total_emitable_copy_num(&copy_nums_true);
    let p = cdbg.is_consistent_copy_num(&copy_nums_true);
    info!("copy-nums total={} consistent={:?}", total, p);
    info!("cycle {}", cdbg.as_cycle_stats());

    for i in 0..cdbg.n_cycles() {
        let is_a = cdbg.is_acceptable(&copy_nums_true, i, true);
        let is_b = cdbg.is_acceptable(&copy_nums_true, i, false);
        let n_rev: u32 = cdbg
            .cycle_components(i)
            .iter()
            .map(|(_, dir)| match dir {
                cycles::CycleDirection::Reverse => 1,
                _ => 0,
            })
            .sum();
        let new_a = cdbg.update_by_cycle(&copy_nums_true, i, true);
        let new_b = cdbg.update_by_cycle(&copy_nums_true, i, false);
        info!(
            "cycle#{} up={} down={} up_ok={} down_ok={} n_rev={}",
            i,
            is_a,
            is_b,
            cdbg.is_consistent_copy_num(&new_a),
            cdbg.is_consistent_copy_num(&new_b),
            n_rev,
        );
    }
}

pub fn sample(dbg_fa: String, length: u32, n_reads: u32, k: usize, seed: u64, param: PHMMParams) {
    let seqs = io::fasta::parse_seqs(&dbg_fa);
    let d = hmm::dbg::DbgPHMM::from_seqs(seqs, k);
    info!("{}", d.dbg.as_degree_stats());
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

pub fn forward(dbg_fa: String, reads_fa: String, k: usize, param: PHMMParams) {
    let seqs = io::fasta::parse_seqs(&dbg_fa);
    let (cdbg, copy_nums) = compressed_dbg::CompressedDBG::from_seqs(seqs, k);
    let phmm = hmm::cdbg::CDbgPHMM::new(&cdbg, copy_nums);

    let reads = io::fasta::parse_seqs(&reads_fa);
    let mut ps: Vec<prob::Prob> = Vec::new();
    for (i, read) in reads.iter().enumerate() {
        let p = phmm.forward_prob(&param, read);
        println!("{}\t{}", i, p.to_log_value());
        ps.push(p);
        // let p = phmm.backward_prob(&param, read);
        // println!("backward prob : {}", p);
    }
    let p_total: prob::Prob = ps.iter().product();
    println!("#total\t{}", p_total.to_log_value());
}

/// Experiments of optimizer
/// 1. construct dbg from reads
/// 2. determine (true) copy_nums from fa
/// 3. optimize
pub fn optimize_with_answer(
    dbg_fa: String,
    reads_fa: String,
    k: usize,
    param: PHMMParams,
    init_temp: f64,
    cooling_rate: f64,
    ave_size: u32,
    std_size: u32,
) {
    let reads = io::fasta::parse_seqs(&reads_fa);
    let (cdbg, _) = compressed_dbg::CompressedDBG::from_seqs(reads, k);

    // let seqs = io::fasta::parse_seqs(&dbg_fa);
}

pub fn optimize(reads_fa: String, k: usize, param: PHMMParams) {
    println!("not implemented!");
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
    compressed_dbg::test();
    // optimizer::cdbg::test();
    // hmm::cdbg::test();
    // distribution::test();
}
