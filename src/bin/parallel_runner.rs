use dbgphmm::e2e::{generate_experiment, ReadType};
use dbgphmm::em::compression::v1::compression_with_depths;
use dbgphmm::em::compression::v3;
use dbgphmm::em::e2e::compression::write_compression_logs;
use dbgphmm::genome;
use dbgphmm::hmmv2::params::PHMMParams;
use itertools::iproduct;
use rayon::prelude::*;
use std::env;
use std::fs::File;
use std::io::Write;
use std::path::Path;

fn run(output_dir: &Path) {
    let mut summary_file = File::create(output_dir.join("summary.txt")).unwrap();
    // writeln!(summary_file, "hoge");

    for (
        (unit_size, n_unit),
        div_init,
        div_hap,
        coverage,
        lambda,
        zero_penalty,
        error_rate,
        k_init,
    ) in iproduct!(
        [(20, 20)],
        [0.1],
        [0.01],
        [10],
        [0.001], // [0.0, 0.0001, 0.001, 0.01, 0.1],
        [-10.0], // [-100.0, -10.0, -5.0, -1.0],
        [0.01],
        [8, 12, 16]
    ) {
        for seed in 0..3 {
            let header = format!(
                "u{unit_size}n{n_unit}di{div_init}dh{div_hap}c{coverage}l{lambda}z{zero_penalty}e{error_rate}seed{seed}ki{k_init}"
            );
            eprintln!("running {}", header);
            let (genome, genome_size) = genome::tandem_repeat_diploid(
                unit_size, n_unit, div_init, seed, seed, div_hap, seed,
            );
            let dataset = generate_experiment(
                genome,
                genome_size,
                seed,
                PHMMParams::uniform(error_rate),
                coverage,
                2000,
                ReadType::FullLengthForHaploid,
                k_init,
                32,
            );

            // v3
            let (new_dbg_v3, logs_v3) = v3::compression(
                &dataset.dbg_raw,
                dataset.reads(),
                &dataset.phmm_params,
                dataset.genome_size(),
                lambda,
                zero_penalty,
                50,
                50,
            );
            let mut task_file =
                File::create(output_dir.join(format!("{}_v3.txt", header))).unwrap();
            write_compression_logs(&mut task_file, &logs_v3, &dataset, &"");

            // v1
            let (new_dbg_v1, logs_v1) = compression_with_depths(
                &dataset.dbg_raw,
                dataset.reads(),
                &dataset.phmm_params,
                &[1.0, 1.0, 2.0, 2.0, 4.0, 4.0],
            );
            let mut task_file =
                File::create(output_dir.join(format!("{}_v1.txt", header))).unwrap();
            for (iteration, log) in logs_v1.iter().enumerate() {
                writeln!(
                    task_file,
                    "{}{}\t{}",
                    header,
                    iteration,
                    log.to_benchmark_string(&dataset)
                );
            }

            let b0 = dataset
                .dbg_true_init
                .benchmark_compression(&dataset, lambda);
            writeln!(summary_file, "{}\tTRUE\t{}", header, b0);
            let br = dataset.dbg_raw.benchmark_compression(&dataset, lambda);
            writeln!(summary_file, "{}\tRAW\t{}", header, br);
            let b3 = new_dbg_v3.benchmark_compression(&dataset, lambda);
            writeln!(summary_file, "{}\tV3\t{}", header, b3);
            let b1 = new_dbg_v1.benchmark_compression(&dataset, lambda);
            writeln!(summary_file, "{}\tV1\t{}", header, b1);
        }
    }
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
