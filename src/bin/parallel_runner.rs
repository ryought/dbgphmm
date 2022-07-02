use dbgphmm::e2e::{generate_dataset, ReadType};
use dbgphmm::em::compression::v3;
use dbgphmm::genome;
use dbgphmm::hmmv2::params::PHMMParams;
use itertools::iproduct;
use rayon::prelude::*;
use std::env;
use std::fs::File;
use std::io::Write;
use std::path::Path;

fn run(output_dir: &Path) {
    let mut summary_file = File::create(output_dir.join("hoge.txt")).unwrap();
    writeln!(summary_file, "hoge");

    /*
    for seed in 0..10 {
        let (genome, genome_size) =
            genome::tandem_repeat_diploid(20, 20, 0.1, seed, seed, 0.01, seed);
        let dataset = generate_dataset(
            genome,
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
        let (new_dbg, logs) = v3::compression(
            &dataset.dbg_raw,
            &dataset.reads,
            &dataset.phmm_params,
            dataset.genome_size,
            lambda,
            50,
            50,
        );
        // inspect_compression_logs(&logs, &genome);
    }
    */
}

fn main() {
    let output_dir = env::args().nth(1);
    match output_dir {
        Some(output_dir) => {
            let output_dir = Path::new(&output_dir);
            if output_dir.exists() {
                // run
                eprintln!("running in {}", output_dir.display());
                run(&output_dir);
            } else {
                eprintln!("ERROR: output dir does not exist.");
            }
        }
        None => {
            eprintln!("ERROR: no output dir was given.");
        }
    }
}
