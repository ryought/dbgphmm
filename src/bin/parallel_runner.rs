use dbgphmm::e2e::{generate_dataset, ReadType};
use dbgphmm::em::compression::v3::compression;
use dbgphmm::genome;
use dbgphmm::hmmv2::params::PHMMParams;
use itertools::iproduct;
use rayon::prelude::*;

fn main() {
    for seed in 0..10 {
        let (genome, genome_size) =
            genome::tandem_repeat_diploid(20, 20, 0.1, seed, seed, 0.01, seed);
        let dataset = generate_dataset(
            &genome,
            genome_size,
            seed,
            PHMMParams::default(),
            10, // coverage
            2000,
            ReadType::FullLength,
            8,
            32,
        );

        let lambda = 0.01;
        let (new_dbg, logs) = compression(
            &dataset.dbg_raw,
            &dataset.reads,
            &dataset.phmm_params,
            genome_size,
            lambda,
            50,
            50,
        );
        // inspect_compression_logs(&logs, &genome);
    }
}
