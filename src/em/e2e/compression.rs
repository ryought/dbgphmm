//!
//! Compression (v3) test
//!
use crate::common::{CopyNum, Genome};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode};
use crate::e2e::{generate_dataset, Experiment, ReadType};
use crate::em::compression::v3::{compression, compression_step, CompressionV3Log};
use std::io::Write;

///
/// utility function for debugging compressionv3 algorithm
///
/// ## TODO
///
/// this should be converted into callback of compression_v3
///
pub fn inspect_compression_logs<N: DbgNode, E: DbgEdge>(
    logs: &[CompressionV3Log<N, E>],
    dataset: &Experiment,
) {
    for (iteration, log) in logs.iter().enumerate() {
        println!("{}\t{}", iteration, log.to_benchmark_string(dataset));
    }
}

pub fn write_compression_logs<N: DbgNode, E: DbgEdge, F: Write>(
    f: &mut F,
    logs: &[CompressionV3Log<N, E>],
    dataset: &Experiment,
    header: &str,
) {
    for (iteration, log) in logs.iter().enumerate() {
        writeln!(
            f,
            "{}{}\t{}",
            header,
            iteration,
            log.to_benchmark_string(dataset)
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
        let lambda = 0.01;
        let (new_dbg, logs) = compression(
            &dataset.dbg_raw,
            &dataset.reads,
            &dataset.phmm_params,
            genome_size,
            lambda,
            -10.0,
            100,
            100,
        );
        inspect_compression_logs(&logs, &dataset);
        println!("dbg_opt={}", new_dbg);
        println!("dbg_tur={}", dataset.dbg_true_init);

        let b0 = dataset
            .dbg_true_init
            .benchmark_compression(&dataset, lambda);
        println!("bench_result={}", b0);
        let br = dataset.dbg_raw.benchmark_compression(&dataset, lambda);
        println!("bench_result={}", br);
        let b = new_dbg.benchmark_compression(&dataset, lambda);
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
            8,
            32,
        );
        println!("genome: {}", genome[0].to_str());
        dataset.show_reads();
        println!("dbg_raw:{}", dataset.dbg_raw);
        println!("dbg_true_init:{}", dataset.dbg_true_init);

        // optimize
        let lambda = 0.0005;
        let (new_dbg, logs) = compression(
            &dataset.dbg_raw,
            &dataset.reads,
            &dataset.phmm_params,
            genome_size,
            lambda,
            -10.0, // clamp
            100,
            100,
        );
        inspect_compression_logs(&logs, &dataset);
        println!("dbg_opt={}", new_dbg);
        println!("dbg_tur={}", dataset.dbg_true_init);

        let b0 = dataset
            .dbg_true_init
            .benchmark_compression(&dataset, lambda);
        println!("bench_result={}", b0);
        let br = dataset.dbg_raw.benchmark_compression(&dataset, lambda);
        println!("bench_result={}", br);
        let b = new_dbg.benchmark_compression(&dataset, lambda);
        println!("bench_result={}", b);
    }

    #[test]
    fn e2e_compression_tandem_repeat_diploid_start_from_true() {
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
            8,
            32,
        );

        // optimize
        let lambda = 0.0005;
        let (new_dbg, logs) = compression(
            &dataset.dbg_true_init,
            &dataset.reads,
            &dataset.phmm_params,
            genome_size,
            lambda,
            -10.0, // clamp
            100,
            100,
        );
        inspect_compression_logs(&logs, &dataset);
        println!("dbg_opt={}", new_dbg);
        println!("dbg_tur={}", dataset.dbg_true_init);

        let b0 = dataset
            .dbg_true_init
            .benchmark_compression(&dataset, lambda);
        println!("bench_result={}", b0);
        let b = new_dbg.benchmark_compression(&dataset, lambda);
        println!("bench_result={}", b);
    }
}
