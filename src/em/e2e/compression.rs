//!
//! Compression (v3) test
//!
use crate::e2e::{generate_dataset, Dataset, ReadType};
use crate::em::compression::v3::compression_step;

///
/// run compression for given dataset
///
pub fn benchmark_compression_v3(dataset: &Dataset) {
    // run compression
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
        let (new_dbg, is_updated, log) = compression_step(
            &dataset.dbg_raw,
            &dataset.reads,
            &dataset.phmm_params,
            genome_size,
            0.0,
        );
        println!("log={}", log);
        println!("dbg_opt={}", new_dbg);
    }
}
