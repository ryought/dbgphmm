//!
//! Tandem repeat haploid/diploid assembly with full-length reads (coverage 10x).
//!

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{sequence_to_string, Genome, Reads, Seq, Sequence};
    use crate::dbg::{Dbg, HashDbg, SimpleDbg};
    use crate::em::compression::{compression, compression_step, compression_with_depths};
    use crate::em::e2e::genome::{generate_tandem_repeat_diploid, generate_tandem_repeat_haploid};
    use crate::em::e2e::runner::generate_full_length_reads_and_dbgs;
    use crate::em::infer;
    use crate::em::scheduler::SchedulerType1;
    use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
    use crate::graph::seq_graph::SeqGraph;
    use crate::hmmv2::params::PHMMParams;
    use crate::hmmv2::sample::{ReadAmount, SampleProfile, StartPoints};
    use crate::kmer::VecKmer;

    fn generate_dataset_haploid(
        unit_size: usize,
        n_unit: usize,
        divergence_init: f64,
        unit_seed: u64,
        hap_seed: u64,
        read_seed: u64,
        phmm_params: PHMMParams,
        coverage: usize,
    ) -> (
        Genome,
        Reads,
        PHMMParams,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
    ) {
        // (1) generate genome and reads
        println!("generating genome");
        let (genome, genome_size) =
            generate_tandem_repeat_haploid(unit_size, n_unit, divergence_init, unit_seed, hap_seed);
        println!("genome hap_a: {}", sequence_to_string(&genome[0]));

        // (2) reads and dbgs
        let (reads, phmm_params, dbg_raw, dbg_true_init, dbg_true) =
            generate_full_length_reads_and_dbgs(
                &genome,
                genome_size,
                read_seed,
                phmm_params,
                coverage,
            );

        (genome, reads, phmm_params, dbg_raw, dbg_true_init, dbg_true)
    }

    fn generate_dataset_diploid(
        unit_seed: u64,
        hap_seed: u64,
        div_seed: u64,
        read_seed: u64,
    ) -> (
        Genome,
        Reads,
        PHMMParams,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
    ) {
        // (1) generate genome and reads
        println!("generating genome");
        let (genome, genome_size) =
            generate_tandem_repeat_diploid(20, 20, 0.1, unit_seed, hap_seed, 0.01, div_seed);
        println!("genome hap_a: {}", sequence_to_string(&genome[0]));
        println!("genome hap_b: {}", sequence_to_string(&genome[1]));

        // (2) reads and dbgs
        let (reads, phmm_params, dbg_raw, dbg_true_init, dbg_true) =
            generate_full_length_reads_and_dbgs(
                &genome,
                genome_size,
                read_seed,
                PHMMParams::default(),
                10,
            );

        (genome, reads, phmm_params, dbg_raw, dbg_true_init, dbg_true)
    }

    #[ignore]
    #[test]
    fn e2e_tandem_repeat_haploid() {
        let (genome, reads, phmm_params, dbg_raw, dbg_true_init, dbg_true) =
            generate_dataset_haploid(20, 20, 0.1, 10, 0, 111, PHMMParams::mid_error(), 30);

        let scheduler = SchedulerType1::new(dbg_raw.k(), dbg_true.k(), 10.0);
        let dbg_infer = infer(&dbg_raw, &reads, &phmm_params, &scheduler, 5);

        println!("{} {}", dbg_true_init.n_traverse_choices(), dbg_true_init);
        println!("{} {}", dbg_true.n_traverse_choices(), dbg_true);
        println!("{} {}", dbg_infer.n_traverse_choices(), dbg_infer);

        let r = dbg_infer.compare(&dbg_true);
        println!("{:?}", r);
        // assert_eq!(r.n_true, 408);
        // assert_eq!(r.n_error, 89);

        let r = dbg_infer.compare_with_seq(&dbg_true, &genome[0]);
        println!("{}", r);

        let p_infer = dbg_infer
            .to_phmm(PHMMParams::mid_error())
            .to_full_prob_parallel(&reads);
        println!("p_infer={}", p_infer);

        let p_true = dbg_true
            .to_phmm(PHMMParams::mid_error())
            .to_full_prob_parallel(&reads);
        println!("p_true={}", p_true);
    }

    #[ignore]
    #[test]
    fn e2e_tandem_repeat_diploid() {
        let (genome, reads, phmm_params, dbg_raw, dbg_true_init, dbg_true) =
            generate_dataset_diploid(10, 0, 2, 111);

        let scheduler = SchedulerType1::new(dbg_raw.k(), dbg_true.k(), 10.0);
        let dbg_infer = infer(&dbg_raw, &reads, &phmm_params, &scheduler, 5);

        println!("{} {}", dbg_true_init.n_traverse_choices(), dbg_true_init);
        println!("{} {}", dbg_true.n_traverse_choices(), dbg_true);
        println!("{} {}", dbg_infer.n_traverse_choices(), dbg_infer);

        let r = dbg_infer.compare(&dbg_true);
        println!("{:?}", r);
        assert_eq!(r.n_true, 732);
        assert_eq!(r.n_error, 0);
    }
}
