use crate::dbg::DBG;
use crate::hmm::base::PHMM;
use crate::hmm::params::PHMMParams;
use crate::hmm::sampler::PHMMSampler;
use crate::kmer::kmer::Kmer;
use crate::*;
use log::{info, warn};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

pub fn generate(length: usize, seed: u64) {
    let v = random_seq::generate(length, seed);
    println!(">randseq");
    println!("{}", std::str::from_utf8(&v).unwrap());
}

pub fn stat(dbg_fa: String, k: usize) {
    let seqs = io::fasta::parse_seqs(&dbg_fa);
    let d = dbg::DbgHash::from_seqs(&seqs, k);
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

pub fn readstat(dbg_fa: String, reads_fa: String, k: usize) {
    let seqs = io::fasta::parse_seqs(&dbg_fa);
    let dbg = dbg::DbgHash::from_seqs(&seqs, k);

    let reads = io::fasta::parse_seqs(&reads_fa);
    let read_dbg = dbg::DbgHash::from_seqs(&reads, k);

    let result = dbg.compare_dbg(&read_dbg);
    println!("{:?}", result);
}

pub fn sample(
    dbg_fa: String,
    length: u32,
    n_reads: u32,
    k: usize,
    seed: u64,
    param: PHMMParams,
    start_from_head: bool,
) {
    let seqs = io::fasta::parse_seqs(&dbg_fa);
    let (cdbg, copy_nums) = compressed_dbg::CompressedDBG::from_seqs(&seqs, k);
    let phmm = hmm::cdbg::CDbgPHMM::new(&cdbg, copy_nums);
    info!("{}", cdbg.as_degree_stats());
    let from = if start_from_head {
        let head = cdbg
            .heads()
            .first()
            .unwrap_or_else(|| panic!("Cannot find head node"))
            .clone();
        Some(head)
    } else {
        None
    };
    info!("from={:?}", from);

    // println!("{}", d.as_dot());
    // println!("{}", d.dbg.as_dot());
    // let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut infos: Vec<hmm::sampler::SampleInfo> = Vec::new();
    for i in 0..n_reads {
        // let seed_for_a_read: u64 = rng.gen();
        let seed_for_a_read = seed + i as u64;
        let (seq, info) = phmm.sample(&param, length, seed_for_a_read, from);
        // output fasta
        let id = format!("{},{}", length, seed_for_a_read);
        io::fasta::dump_seq(&id, &seq, Some(&info.to_string()));
        // store info in vec
        infos.push(info);
    }
    info!("{:?}", hmm::sampler::sum_sample_infos(&infos));
}

pub fn forward(dbg_fa: String, reads_fa: String, k: usize, param: PHMMParams, parallel: bool) {
    let seqs = io::fasta::parse_seqs(&dbg_fa);
    let (cdbg, copy_nums) = compressed_dbg::CompressedDBG::from_seqs(&seqs, k);
    let phmm = hmm::cdbg::CDbgPHMM::new(&cdbg, copy_nums);

    let reads = io::fasta::parse_seqs(&reads_fa);

    if parallel {
        let p_total: prob::Prob = reads
            .par_iter()
            .map(|read| phmm.forward_prob(&param, read))
            .product();
        println!("#total\t{}", p_total.to_log_value());
    } else {
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
    n_iteration: u64,
    ave_size: u32,
    std_size: u32,
    prior_only: bool,
    start_from_true_copy_nums: bool,
    parallel: bool,
) {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
    let a = optimizer::base::Annealer::new(init_temp, cooling_rate);

    let reads = io::fasta::parse_seqs(&reads_fa);
    let (cdbg, _) = compressed_dbg::CompressedDBG::from_seqs(&reads, k);

    let seqs = io::fasta::parse_seqs(&dbg_fa);
    let copy_nums_true = cdbg
        .true_copy_nums_from_seqs(&seqs, k)
        .unwrap_or_else(|| panic!("True copy_nums is not in read cdbg"));
    let true_size = cdbg.total_emitable_copy_num(&copy_nums_true);
    info!("true_size={}", true_size);

    let init_state = optimizer::cdbg::CDbgState::init(
        &cdbg,
        true_size,
        std_size,
        if prior_only { None } else { Some(&reads) },
        param.clone(),
        parallel,
    );
    let cycle_vec_true = cdbg.cycle_vec_from_copy_nums(&copy_nums_true);
    let true_state = optimizer::cdbg::CDbgState::new(
        &cdbg,
        copy_nums_true.clone(),
        cycle_vec_true,
        true_size,
        std_size,
        if prior_only { None } else { Some(&reads) },
        param.clone(),
        parallel,
    );

    if start_from_true_copy_nums {
        let history = a.run_with_log(&mut rng, true_state, n_iteration);
        let copy_nums_final = &history.last().unwrap().copy_nums;
        for (i, seq) in cdbg.to_seqs(copy_nums_final).iter().enumerate() {
            let id = format!("{}", i);
            io::fasta::dump_seq(&id, &seq, None);
        }
    } else {
        let history = a.run_with_log(&mut rng, init_state, n_iteration);
        a.run_with_log(&mut rng, true_state, 1);
        // println!("{:?}", history.last().unwrap().copy_nums);
        // println!("{:?}", copy_nums_true);
        let copy_nums_final = &history.last().unwrap().copy_nums;
        for (i, seq) in cdbg.to_seqs(copy_nums_final).iter().enumerate() {
            let id = format!("{}", i);
            io::fasta::dump_seq(&id, &seq, None);
        }
    }

    for (i, seq) in cdbg.to_seqs(&copy_nums_true).iter().enumerate() {
        let id = format!("t{}", i);
        io::fasta::dump_seq(&id, &seq, None);
    }
}

pub fn optimize(reads_fa: String, k: usize, param: PHMMParams) {
    println!("not implemented!");
}

pub fn sandbox() {
    let mut reads: Vec<Vec<u8>> = Vec::new();
    // reads.push(b"ATCGATTCGATCGATTCGATAGATCG".to_vec());
    reads.push(b"AGGCTAGTAAAAAAAAAAAAAATCGATCTTTCGATCG".to_vec());
    reads.push(b"GGATAGTTCGATCTG".to_vec());
    reads.push(b"GGCTAGTTCGATCGG".to_vec());
    let (cdbg, copy_nums) = compressed_dbg::CompressedDBG::from_seqs(&reads, 8);
    // println!("{}", cdbg.as_dot_with_copy_nums(&copy_nums));
    for (i, seq) in cdbg.to_seqs(&copy_nums).iter().enumerate() {
        let id = format!("{}", i);
        io::fasta::dump_seq(&id, &seq, None);
    }
}
