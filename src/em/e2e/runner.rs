//!
//! read and dbg generation functions
//!
use crate::common::{sequence_to_string, Genome, Reads, Seq, Sequence};
use crate::dbg::compare::{CompareResult, CompareWithSeqResult};
use crate::dbg::{Dbg, HashDbg, SimpleDbg};
use crate::em::infer;
use crate::em::scheduler::SchedulerType1;
use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::sample::{ReadAmount, SampleProfile, StartPoints};
use crate::kmer::VecKmer;

pub fn benchmark_em_steps(
    dbg_raw: &SimpleDbg<VecKmer>,
    dbg_true: &SimpleDbg<VecKmer>,
    reads: &Reads,
    genome: &Genome,
    phmm_params: &PHMMParams,
    coverage: f64,
) {
    let scheduler = SchedulerType1::new(dbg_raw.k(), dbg_true.k(), coverage);
    let (dbg_infer, logs) = infer(dbg_raw, reads, phmm_params, &scheduler, 5);
    for (i, log) in logs.iter().enumerate() {
        println!("iter#{}\t{:?}", i, log);
    }
}

pub fn benchmark(
    dbg_raw: &SimpleDbg<VecKmer>,
    dbg_true: &SimpleDbg<VecKmer>,
    reads: &Reads,
    genome: &Genome,
    phmm_params: &PHMMParams,
    coverage: f64,
) -> (
    SimpleDbg<VecKmer>,
    CompareResult,
    Vec<CompareWithSeqResult<VecKmer>>,
) {
    let scheduler = SchedulerType1::new(dbg_raw.k(), dbg_true.k(), coverage);
    let (dbg_infer, _) = infer(dbg_raw, reads, phmm_params, &scheduler, 5);

    println!("dbg_raw=\n{}", dbg_raw);
    println!("{}", dbg_raw.n_traverse_choices());
    println!("dbg_infer=\n{}", dbg_infer);
    println!("{}", dbg_infer.n_traverse_choices());
    println!("dbg_true=\n{}", dbg_true);
    println!("{}", dbg_true.n_traverse_choices());
    for (i, haplotype) in genome.iter().enumerate() {
        println!("genome{}=\n{}", i, sequence_to_string(haplotype));
    }

    let r = dbg_infer.compare(&dbg_true);
    println!("{:?}", r);

    let p_infer = dbg_infer
        .to_phmm(phmm_params.clone())
        .to_full_prob_parallel(reads);
    println!("p_infer={}", p_infer);

    let p_true = dbg_true
        .to_phmm(phmm_params.clone())
        .to_full_prob_parallel(reads);
    println!("p_true={}", p_true);

    let mut v = Vec::new();
    for hap in genome {
        let rs = dbg_infer.compare_with_seq(&dbg_true, hap);
        println!("{}", rs);
        v.push(rs);
    }

    (dbg_infer, r, v)
}

pub enum ReadType {
    FullLength,
    Fragment,
}

pub fn generate_reads_and_dbgs(
    genome: &Genome,
    genome_size: usize,
    read_seed: u64,
    phmm_params: PHMMParams,
    coverage: usize,
    read_length: usize,
    read_type: ReadType,
    k_init: usize,
    k_target: usize,
) -> (
    Reads,
    PHMMParams,
    SimpleDbg<VecKmer>,
    SimpleDbg<VecKmer>,
    SimpleDbg<VecKmer>,
) {
    let g = GenomeGraph::from_seqs(genome);
    let profile = match read_type {
        ReadType::Fragment => ReadProfile {
            has_revcomp: true,
            sample_profile: SampleProfile {
                read_amount: ReadAmount::TotalBases(genome_size * coverage),
                seed: read_seed,
                length: read_length,
                start_points: StartPoints::Random,
                endable: false,
            },
            phmm_params: phmm_params.clone(),
        },
        ReadType::FullLength => ReadProfile {
            has_revcomp: false,
            sample_profile: SampleProfile {
                read_amount: ReadAmount::Count(coverage),
                seed: read_seed,
                length: read_length,
                start_points: StartPoints::AllStartPoints,
                endable: false,
            },
            phmm_params: phmm_params.clone(),
        },
    };
    let pos_reads = g.sample_positioned_reads(&profile);
    g.show_coverage(&pos_reads);
    let reads = pos_reads.to_reads(true);
    for (i, read) in reads.iter().enumerate() {
        println!("read#{}\t{}", i, read.to_str());
    }

    let dbg_raw: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k_init, &reads);

    // (4) compare with true dbg
    let dbg_true_init: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k_init, genome);

    // (5) true k=50 (read length)
    let dbg_true: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k_target, genome);

    (reads, phmm_params, dbg_raw, dbg_true_init, dbg_true)
}

pub fn generate_full_length_reads_and_dbgs(
    genome: &Genome,
    genome_size: usize,
    read_seed: u64,
    phmm_params: PHMMParams,
    coverage: usize,
) -> (
    Reads,
    PHMMParams,
    SimpleDbg<VecKmer>,
    SimpleDbg<VecKmer>,
    SimpleDbg<VecKmer>,
) {
    let g = GenomeGraph::from_seqs(genome);
    let profile = ReadProfile {
        has_revcomp: true,
        sample_profile: SampleProfile {
            read_amount: ReadAmount::TotalBases(genome_size * coverage),
            seed: read_seed,
            length: 1000,
            start_points: StartPoints::AllStartPoints,
            endable: false,
        },
        phmm_params: phmm_params.clone(),
    };
    let pos_reads = g.sample_positioned_reads(&profile);
    g.show_coverage(&pos_reads);
    let reads = pos_reads.to_reads(true);
    for read in reads.iter() {
        println!("{}", read.to_str());
    }

    let k: usize = 8;
    let dbg_raw: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k, &reads);

    // (4) compare with true dbg
    let dbg_true_init: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k, genome);

    // (5) true k=50 (read length)
    let dbg_true: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(100, genome);

    (reads, phmm_params, dbg_raw, dbg_true_init, dbg_true)
}
