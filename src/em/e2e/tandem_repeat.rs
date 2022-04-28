//!
//! Tandem repeat haploid/diploid assembly with full-length reads (coverage 10x).
//!
//! ## Ignored tests
//!
//! * e2e_tandem_repeat_haploid_error_2
//! * e2e_tandem_repeat_haploid_error_1
//! * e2e_tandem_repeat_diploid
//!

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{sequence_to_string, Genome, Reads, Seq, Sequence};
    use crate::dbg::{Dbg, HashDbg, SimpleDbg};
    use crate::em::compression::{compression, compression_step, compression_with_depths};
    use crate::em::e2e::genome::{generate_tandem_repeat_diploid, generate_tandem_repeat_haploid};
    use crate::em::e2e::runner::{benchmark, generate_full_length_reads_and_dbgs};
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
    fn e2e_tandem_repeat_haploid_error_2() {
        let (genome, reads, phmm_params, dbg_raw, dbg_true_init, dbg_true) =
            generate_dataset_haploid(20, 20, 0.1, 10, 0, 111, PHMMParams::mid_error_2(), 10);
        let (_, r, s) = benchmark(&dbg_raw, &dbg_true, &reads, &genome, &phmm_params, 10.0);
        assert_eq!(r.n_true, 248);
        assert_eq!(r.n_error, 249);
        for (i, result) in s[0].0.iter().enumerate() {
            if i <= 221 {
                assert!(result.is_correct());
            } else if i <= 470 {
                assert!(!result.is_correct());
            } else {
                assert!(result.is_correct());
            }
        }
    }

    #[ignore]
    #[test]
    fn e2e_tandem_repeat_haploid_error_1() {
        let (genome, reads, phmm_params, dbg_raw, dbg_true_init, dbg_true) =
            generate_dataset_haploid(20, 20, 0.1, 10, 0, 111, PHMMParams::default(), 10);
        let (_, r, s) = benchmark(&dbg_raw, &dbg_true, &reads, &genome, &phmm_params, 10.0);
        assert_eq!(r.n_true, 408);
        assert_eq!(r.n_error, 89);
        for (i, result) in s[0].0.iter().enumerate() {
            if i <= 397 {
                assert!(result.is_correct());
            }
        }
    }

    #[ignore]
    #[test]
    fn e2e_tandem_repeat_haploid_error_v2_1() {
        let (genome, reads, phmm_params, dbg_raw, dbg_true_init, dbg_true) =
            generate_dataset_haploid(20, 20, 0.1, 10, 0, 111, PHMMParams::default(), 15);
        let (_, r, s) = benchmark(&dbg_raw, &dbg_true, &reads, &genome, &phmm_params, 15.0);
        assert_eq!(r.n_true, 362);
        assert_eq!(r.n_error, 135);
    }

    #[ignore]
    #[test]
    fn e2e_tandem_repeat_diploid() {
        let (genome, reads, phmm_params, dbg_raw, dbg_true_init, dbg_true) =
            generate_dataset_diploid(10, 0, 2, 111);
        let (_, r, _) = benchmark(&dbg_raw, &dbg_true, &reads, &genome, &phmm_params, 10.0);
        assert_eq!(r.n_true, 732);
        assert_eq!(r.n_error, 0);
    }
}
