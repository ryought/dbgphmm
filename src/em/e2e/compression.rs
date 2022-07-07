//!
//! Compression (v3) test
//!
use crate::common::{CopyNum, Genome};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode};
use crate::e2e::{generate_dataset, Dataset, ReadType};
use crate::em::compression::v3::{compression, compression_step, CompressionV3Log};
use std::io::Write;

pub fn inspect_compression_logs<N: DbgNode, E: DbgEdge>(
    logs: &[CompressionV3Log<N, E>],
    genome: &Genome,
) {
    write_compression_logs(&mut std::io::stdout(), logs, genome, &"inspect\t");
}

pub fn write_compression_logs<N: DbgNode, E: DbgEdge, F: Write>(
    f: &mut F,
    logs: &[CompressionV3Log<N, E>],
    genome: &Genome,
    header: &str,
) {
    for (iteration, log) in logs.iter().enumerate() {
        // println!("log={}", log);
        let kh = log.dbg.kmer_hists_from_seqs(genome);
        let ke = log.dbg.check_kmer_existence_with_seqs(genome);
        writeln!(
            f,
            "{}{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            header,
            iteration,
            log.p.to_log_value(),
            log.q0,
            log.q1,
            log.cost_diff,
            log.dbg.genome_size(),
            kh.n_missed_kmers(),
            kh.n_under_estimated_kmers(),
            ke,
            kh,
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
    use crate::min_flow::utils::DEFAULT_CLAMP_VALUE;

    #[test]
    fn e2e_compression_simple() {
        // data generation
        let (genome, genome_size) = genome::simple(200, 5);
        let dataset = generate_dataset(
            genome.clone(),
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
            0.1,
            DEFAULT_CLAMP_VALUE,
            10,
            10,
        );
        inspect_compression_logs(&logs, &genome);
        println!("dbg_opt={}", new_dbg);
        println!("dbg_tur={}", dataset.dbg_true_init);

        let b0 = dataset.dbg_true_init.benchmark_compression(&dataset);
        println!("bench_result={}", b0);
        let br = dataset.dbg_raw.benchmark_compression(&dataset);
        println!("bench_result={}", br);
        let b = new_dbg.benchmark_compression(&dataset);
        println!("bench_result={}", b);
    }

    #[test]
    fn e2e_compression_tandem_repeat_diploid() {
        // data generation
        let (genome, genome_size) = genome::tandem_repeat_diploid(20, 20, 0.1, 0, 0, 0.01, 0);
        let dataset = generate_dataset(
            genome.clone(),
            genome_size,
            0,
            PHMMParams::default(),
            10, // coverage
            2000,
            ReadType::FullLength,
            16,
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
            0.0001, // lambda
            -10.0,  // clamp
            50,
            50,
        );
        inspect_compression_logs(&logs, &genome);
        println!("dbg_opt={}", new_dbg);
        println!("dbg_tur={}", dataset.dbg_true_init);

        let b0 = dataset.dbg_true_init.benchmark_compression(&dataset);
        println!("bench_result={}", b0);
        let br = dataset.dbg_raw.benchmark_compression(&dataset);
        println!("bench_result={}", br);
        let b = new_dbg.benchmark_compression(&dataset);
        println!("bench_result={}", b);
    }
}
