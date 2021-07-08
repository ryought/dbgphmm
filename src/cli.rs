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
    info!("{:?}", d.as_degree_stats());
    let cdbg = compressed_dbg::CompressedDBG::from(&d, k);
    info!("{:?}", cdbg.as_dbg_stats());
    info!("{:?}", cdbg.as_degree_stats());
    let copy_nums_true: Vec<u32> = cdbg
        .iter_nodes()
        .map(|v| {
            let kmer = cdbg.kmer(&v);
            d.find(kmer)
        })
        .collect();

    info!("{:?}", cdbg.as_copy_num_stats(&copy_nums_true));
    // let cn_zero: Vec<u32> = vec![0; cdbg.n_kmers()];
    // info!("{:?}", cdbg.as_copy_num_stats(&cn_zero));
    // info!("{:?}", copy_nums_true);
    info!("cycle {}", cdbg.as_cycle_histogram());
    info!("{:?}", cdbg.as_cycle_summary_stats());

    for i in 0..cdbg.n_cycles() {
        info!("{:?}", cdbg.as_cycle_stats(i));
    }

    let all_stats = cdbg.as_all_stats(&copy_nums_true);
    let json = serde_json::to_string_pretty(&all_stats).unwrap();
    // info!("{:#?}", cdbg.as_all_stats(&copy_nums_true));
    println!("{}", json);
}

pub fn compare(self_dbg_fa: String, other_dbg_fa: String, k: usize) {
    let self_seqs = io::fasta::parse_seqs(&self_dbg_fa);
    let self_dbg = dbg::DbgHash::from_seqs(&self_seqs, k);

    let other_seqs = io::fasta::parse_seqs(&other_dbg_fa);
    let other_dbg = dbg::DbgHash::from_seqs(&other_seqs, k);

    let result = self_dbg.compare_dbg(&other_dbg);
    let json = serde_json::to_string_pretty(&result).unwrap();
    println!("{}", json);
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
    info!("{:?}", cdbg.as_degree_stats());
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
    dump_seqs: bool,
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
        if dump_seqs {
            for (i, seq) in cdbg.to_seqs(copy_nums_final).iter().enumerate() {
                let id = format!("{}", i);
                io::fasta::dump_seq(&id, &seq, None);
            }
        }
    } else {
        let history = a.run_with_log(&mut rng, init_state, n_iteration);
        a.run_with_log(&mut rng, true_state, 1);
        // println!("{:?}", history.last().unwrap().copy_nums);
        // println!("{:?}", copy_nums_true);
        let copy_nums_final = &history.last().unwrap().copy_nums;
        if dump_seqs {
            for (i, seq) in cdbg.to_seqs(copy_nums_final).iter().enumerate() {
                let id = format!("{}", i);
                io::fasta::dump_seq(&id, &seq, None);
            }
        }
    }

    if dump_seqs {
        for (i, seq) in cdbg.to_seqs(&copy_nums_true).iter().enumerate() {
            let id = format!("t{}", i);
            io::fasta::dump_seq(&id, &seq, None);
        }
    }
}

pub fn optimize(reads_fa: String, k: usize, param: PHMMParams) {
    println!("not implemented!");
}

pub fn sandbox() {
    stats::test();
    /*
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
    */
}
