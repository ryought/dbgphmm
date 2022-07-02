//!
//! Compression (v3) test
//!
use crate::common::{CopyNum, Genome};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode};
use crate::e2e::{generate_dataset, Dataset, ReadType};
use crate::em::compression::v3::{compression, compression_step, CompressionV3Log};

pub fn inspect_compression_logs<N: DbgNode, E: DbgEdge>(
    logs: &[CompressionV3Log<N, E>],
    genome: &Genome,
) {
    for (iteration, log) in logs.iter().enumerate() {
        // println!("log={}", log);
        let kh = log.dbg.kmer_hists_from_seqs(genome);
        let ke = log.dbg.check_kmer_existence_with_seqs(genome);
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            iteration,
            log.p.to_log_value(),
            log.q0,
            log.q1,
            log.cost_diff,
            log.dbg.genome_size(),
            kh.n_missed_kmers(),
            ke,
            kh,
            // log.dbg,
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
            10,
            10,
        );
        inspect_compression_logs(&logs, &genome);
        println!("dbg_opt={}", new_dbg);
        println!("dbg_tur={}", dataset.dbg_true_init);
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
        let (new_dbg, logs) = compression(
            &dataset.dbg_raw,
            &dataset.reads,
            &dataset.phmm_params,
            genome_size,
            0.01,
            50,
            50,
        );
        inspect_compression_logs(&logs, &genome);
        println!("dbg_opt={}", new_dbg);
        println!("dbg_tur={}", dataset.dbg_true_init);
    }
}
