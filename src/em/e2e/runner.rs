//!
//! read and dbg generation functions
//!
use crate::common::collection::genome_size;
use crate::common::{sequence_to_string, Genome, Reads, Seq, Sequence};
use crate::dbg::compare::{CompareResult, CompareWithSeqResult};
use crate::dbg::dbg::{DbgEdge, DbgNode};
use crate::dbg::{Dbg, HashDbg, SimpleDbg};
use crate::e2e::{Experiment, ReadType};
use crate::em::e2e::compression::write_compression_logs;
use crate::em::infer;
use crate::em::scheduler::SchedulerType1;
use crate::em::TaskLog;
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
pub fn show_logs<N: DbgNode, E: DbgEdge>(task_logs: &[TaskLog<N, E>], dataset: &Experiment) {
    // header
    println!("iter\ttype\tstep\tprob\tmin_flow\tgenome_size\tcompare\tdbg\t");
    // body
    for (iteration, task_log) in task_logs.iter().enumerate() {
        println!(
            "{}",
            task_log.to_benchmark_string_with_header(&dataset, &format!("{}\t", iteration))
        );
    }
}

pub fn benchmark(
    dataset: &Experiment,
    coverage: f64,
) -> (
    SimpleDbg<VecKmer>,
    CompareResult,
    Vec<CompareWithSeqResult<VecKmer>>,
) {
    let scheduler = SchedulerType1::new(dataset.dbg_raw.k(), dataset.dbg_true.k(), coverage);
    let genome_size = dataset.genome_size();
    let (dbg_infer, logs) = infer(
        &dataset.dbg_raw,
        &dataset.reads().clone().to_reads(),
        &dataset.phmm_params,
        &scheduler,
        genome_size,
        5,
    );

    println!("dbg_raw=\n{}", dataset.dbg_raw);
    println!("{}", dataset.dbg_raw.n_traverse_choices());
    println!("dbg_infer=\n{}", dbg_infer);
    println!("{}", dbg_infer.n_traverse_choices());
    println!("dbg_true=\n{}", dataset.dbg_true);
    println!("{}", dataset.dbg_true.n_traverse_choices());
    for (i, haplotype) in dataset.genome().iter().enumerate() {
        println!("genome{}=\n{}", i, sequence_to_string(haplotype));
    }

    let r = dbg_infer.compare(&dataset.dbg_true);
    println!("{:?}", r);

    let p_infer = dbg_infer
        .to_phmm(dataset.phmm_params.clone())
        .to_full_prob_parallel(dataset.reads());
    println!("p_infer={}", p_infer);

    let p_true = dataset
        .dbg_true
        .to_phmm(dataset.phmm_params.clone())
        .to_full_prob_parallel(dataset.reads());
    println!("p_true={}", p_true);

    let mut v = Vec::new();
    for hap in dataset.genome() {
        let rs = dbg_infer.compare_with_seq(&dataset.dbg_true, hap);
        println!("{}", rs);
        v.push(rs);
    }

    show_logs(&logs, &dataset);

    (dbg_infer, r, v)
}
