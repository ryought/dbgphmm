//!
//! Compression (v3) test
//!
use crate::common::Genome;
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode};
use crate::e2e::{generate_dataset, Dataset, ReadType};
use crate::em::compression::v3::{compression, compression_step, CompressionV3Log};

///
/// run compression for given dataset
///
pub fn benchmark_compression_v3(dataset: &Dataset) {
    // run compression
}

pub fn inspect_compression_logs<N: DbgNode, E: DbgEdge>(
    logs: &[CompressionV3Log<N, E>],
    genome: &Genome,
) {
    for (iteration, log) in logs.iter().enumerate() {
        // println!("log={}", log);
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            iteration,
            log.p.to_log_value(),
            log.q0,
            log.q1,
            log.cost_diff,
            log.dbg.genome_size(),
            log.dbg.kmer_hists_from_seqs(genome),
            log.dbg,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Seq;
    use crate::genome;
    use crate::hmmv2::params::PHMMParams;
    #[test]
    fn e2e_compression_simple() {
        // data generation
        let (genome, genome_size) = genome::simple(1000, 5);
        let dataset = generate_dataset(
            &genome,
            genome_size,
            0,
            PHMMParams::default(),
            10, // coverage
            2000,
            ReadType::FullLength,
            8,
            32,
        );
        println!("genome: {}", genome[0].to_str());
        dataset.show_reads();
        println!("dbg_raw:{}", dataset.dbg_raw);
        println!("dbg_true_init:{}", dataset.dbg_true_init);

        // optimize
        let (new_dbg, logs) = compression(
            &dataset.dbg_raw,
            &dataset.reads,
            &dataset.phmm_params,
            genome_size,
            1.0,
            10,
            10,
        );
        inspect_compression_logs(&logs, &genome);
        println!("dbg_opt={}", new_dbg);
        println!("dbg_tur={}", dataset.dbg_true_init);
    }
}
