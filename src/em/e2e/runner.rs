//!
//! read and dbg generation functions
//!
use crate::common::collection::genome_size;
use crate::common::{sequence_to_string, Genome, Reads, Seq, Sequence};
use crate::dbg::compare::{CompareResult, CompareWithSeqResult};
use crate::dbg::dbg::{DbgEdge, DbgNode};
use crate::dbg::{Dbg, HashDbg, SimpleDbg};
use crate::e2e::{Dataset, ReadType};
use crate::em::e2e::compression::write_compression_logs;
use crate::em::infer;
use crate::em::scheduler::{SchedulerType1, TaskLog};
use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::sample::{ReadAmount, SampleProfile, StartPoints};
use crate::kmer::VecKmer;
use std::io::Write;

///
/// WillDeprecate
///
/// show a TaskLog list with true genome
///
pub fn show_logs<N: DbgNode, E: DbgEdge>(task_logs: &[TaskLog<N, E>], genome: &Genome) {
    // header
    println!("iter\ttype\tstep\tprob\tmin_flow\tgenome_size\tcompare\tdbg\t");

    // body
    for (iteration, task_log) in task_logs.iter().enumerate() {
        match task_log {
            TaskLog::Compression(logs) => {
                for (step, log) in logs.iter().enumerate() {
                    println!(
                        "{}\tC\t{}\t{}\t{}\t{}\t{}\t{}",
                        iteration,
                        step,
                        log.full_prob.to_log_value(),
                        log.min_flow_score,
                        log.dbg.genome_size(),
                        log.dbg.kmer_hists_from_seqs(genome),
                        log.dbg,
                    );
                }
            }
            TaskLog::Extension(logs) => {
                for (step, log) in logs.iter().enumerate() {
                    println!(
                        "{}\tE\t{}\t{}\t{}\t{}\t{}\t{}",
                        iteration,
                        step,
                        match log.full_prob {
                            Some(p) => p.to_log_value().to_string(),
                            None => "-".to_string(),
                        },
                        log.min_flow_cost,
                        log.dbg.genome_size(),
                        log.dbg.kmer_hists_from_seqs(genome),
                        log.dbg
                    );
                }
            }
            _ => panic!(),
        }
    }
}

pub fn benchmark_em_steps(
    dbg_raw: &SimpleDbg<VecKmer>,
    dbg_true: &SimpleDbg<VecKmer>,
    reads: &Reads,
    genome: &Genome,
    phmm_params: &PHMMParams,
    coverage: f64,
) {
    let scheduler = SchedulerType1::new(dbg_raw.k(), dbg_true.k(), coverage);
    let genome_size = genome_size(genome);
    let (dbg_infer, logs) = infer(dbg_raw, reads, phmm_params, &scheduler, genome_size, 5);
    show_logs(&logs, genome);
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
    let genome_size = genome_size(genome);
    let (dbg_infer, logs) = infer(dbg_raw, reads, phmm_params, &scheduler, genome_size, 5);

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

    show_logs(&logs, genome);

    (dbg_infer, r, v)
}
