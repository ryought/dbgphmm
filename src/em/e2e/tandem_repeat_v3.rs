//!
//! Tandem repeat haploid/diploid assembly test with v3 compression
//!

#[cfg(test)]
mod tests {
    use crate::common::Seq;
    use crate::e2e::{generate_dataset, Dataset, ReadType};
    use crate::em::infer;
    use crate::em::scheduler::SchedulerType1;
    use crate::genome;
    use crate::hmmv2::params::PHMMParams;
    #[test]
    fn e2e_v3_tandem_repeat_diploid() {
        // data generation
        let (genome, genome_size) = genome::tandem_repeat_diploid(20, 20, 0.1, 0, 0, 0.01, 0);
        let coverage = 10;
        let dataset = generate_dataset(
            genome.clone(),
            genome_size,
            0,
            PHMMParams::default(),
            coverage,
            2000,
            ReadType::FullLength,
            16,
            32,
        );
        println!("genome: {}", genome[0].to_str());
        dataset.show_reads();
        println!("dbg_raw:{}", dataset.dbg_raw);
        println!("dbg_true_init:{}", dataset.dbg_true_init);

        let scheduler = SchedulerType1::new_v3(
            dataset.dbg_raw.k(),
            dataset.dbg_true.k(),
            10,
            -0.001,
            0.1,
            -100.0,
            -100.0,
        );
        let (dbg_infer, logs) = infer(
            &dataset.dbg_raw,
            &dataset.reads,
            &dataset.phmm_params,
            &scheduler,
            genome_size,
            50,
        );
    }
}
