//!
//! Tandem repeat haploid/diploid assembly test with v3 compression
//!

#[cfg(test)]
mod tests {
    use crate::common::Seq;
    use crate::e2e::{generate_dataset, Experiment, ReadType};
    use crate::em::infer_with_on_iteration;
    use crate::em::scheduler::SchedulerType1;
    use crate::genome;
    use crate::hmmv2::params::PHMMParams;
    #[test]
    fn e2e_v1_vs_v3_simple() {
        // data generation
        // let (genome, genome_size) = genome::diploid(200, 3, 0.1, 0);
        let (genome, genome_size) = genome::simple(200, 5);
        let coverage = 10;
        let dataset = generate_dataset(
            genome.clone(),
            genome_size,
            0,
            PHMMParams::uniform(0.01),
            coverage,
            2000,
            ReadType::FullLength,
            8,
            32,
        );
        println!("genome: {}", genome[0].to_str());
        dataset.show_reads();
        println!("dbg_raw:{}", dataset.dbg_raw);
        println!("dbg_true_init:{}", dataset.dbg_true_init);

        /*
        let scheduler =
            SchedulerType1::new(dataset.dbg_raw.k(), dataset.dbg_true.k(), coverage as f64);
        */
        let scheduler = SchedulerType1::new_v3(
            dataset.dbg_raw.k(),
            dataset.dbg_true.k(),
            10,
            0.0001,
            0.01,
            -1000.0,
            -1000.0,
        );

        let (dbg_infer, logs) = infer_with_on_iteration(
            &dataset.dbg_raw,
            &dataset.reads,
            &dataset.phmm_params,
            &scheduler,
            genome_size,
            50,
            |iteration, task, task_log, dbg| {
                println!(
                    "{}",
                    task_log.to_benchmark_string_with_header(&dataset, &format!("{}\t", iteration))
                );
            },
        );
    }
    #[test]
    fn e2e_v1_vs_v3_tandem_repeat_diploid() {
        // data generation
        let (genome, genome_size) = genome::tandem_repeat_diploid(20, 20, 0.1, 0, 0, 0.01, 0);
        let coverage = 10;
        let dataset = generate_dataset(
            genome.clone(),
            genome_size,
            0,
            PHMMParams::uniform(0.003),
            coverage,
            2000,
            ReadType::FullLength,
            8,
            32,
        );
        println!("genome: {}", genome[0].to_str());
        dataset.show_reads();
        println!("dbg_raw:{}", dataset.dbg_raw);
        println!("dbg_true_init:{}", dataset.dbg_true_init);

        /*
        let scheduler = SchedulerType1::new_v3(
            dataset.dbg_raw.k(),
            dataset.dbg_true.k(),
            10,
            0.001,
            0.1,
            -100.0,
            -100.0,
        );
        */
        let scheduler =
            SchedulerType1::new(dataset.dbg_raw.k(), dataset.dbg_true.k(), coverage as f64);

        let (dbg_infer, logs) = infer_with_on_iteration(
            &dataset.dbg_raw,
            &dataset.reads,
            &dataset.phmm_params,
            &scheduler,
            genome_size,
            50,
            |iteration, task, task_log, dbg| {
                println!(
                    "{}",
                    task_log.to_benchmark_string_with_header(&dataset, &format!("{}\t", iteration))
                );
            },
        );
    }
}
