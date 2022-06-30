//!
//! Compression (v3) test
//!
use crate::e2e::{generate_dataset, Dataset, ReadType};

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
        let (genome, genome_size) = genome::simple(200, 5);
        println!("genome: {}", genome[0].to_str());

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

        dataset.show_reads();
        println!("dbg_raw:{}", dataset.dbg_raw);
        println!("dbg_true_init:{}", dataset.dbg_true_init);
    }
}
