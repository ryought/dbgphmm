//!
//! Tandem repeat haploid/diploid assembly with full-length reads (coverage 10x).
//!

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{sequence_to_string, Genome, Reads, Seq, Sequence};
    use crate::dbg::{Dbg, HashDbg, SimpleDbg};
    use crate::em::compression::{compression, compression_step, compression_with_depths};
    use crate::em::infer;
    use crate::em::scheduler::SchedulerType1;
    use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
    use crate::graph::seq_graph::SeqGraph;
    use crate::hmmv2::params::PHMMParams;
    use crate::hmmv2::sample::{ReadAmount, SampleProfile, StartPoints};
    use crate::kmer::VecKmer;
    use crate::random_seq::{generate, random_mutation, tandem_repeat, MutationProfile};

    fn generate_genome_haploid(unit_seed: u64, hap_seed: u64) -> (Genome, usize) {
        let unit_size = 20;
        let n_unit = 20;
        let genome_size = unit_size * n_unit;
        let unit = generate(unit_size, unit_seed);
        let tandem_repeat = tandem_repeat(&unit, n_unit);
        let (hap_a, _) = random_mutation(&tandem_repeat, MutationProfile::uniform(0.1), hap_seed);
        println!("{}", tandem_repeat.to_str());
        println!("{}", hap_a.to_str());
        (vec![hap_a], genome_size)
    }

    fn generate_genome_diploid(unit_seed: u64, hap_seed: u64, div_seed: u64) -> (Genome, usize) {
        let (mut hap, hap_genome_size) = generate_genome_haploid(unit_seed, hap_seed);
        let hap_a = hap.remove(0);
        let (hap_b, _) = random_mutation(&hap_a, MutationProfile::uniform(0.01), div_seed);
        let genome_size = hap_genome_size * 2;
        println!("{}", hap_a.to_str());
        println!("{}", hap_b.to_str());
        (vec![hap_a, hap_b], genome_size)
    }

    fn generate_full_length_reads_and_dbgs(
        genome: &Genome,
        genome_size: usize,
        read_seed: u64,
    ) -> (
        Reads,
        PHMMParams,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
    ) {
        println!("generating reads");
        let g = GenomeGraph::from_seqs(genome);
        let phmm_params = PHMMParams::default();
        let profile = ReadProfile {
            has_revcomp: true,
            sample_profile: SampleProfile {
                read_amount: ReadAmount::TotalBases(genome_size * 10),
                seed: read_seed,
                length: 1000,
                start_points: StartPoints::AllStartPoints,
                endable: false,
            },
            phmm_params: phmm_params.clone(),
        };
        let pos_reads = g.sample_positioned_reads(&profile);
        g.show_coverage(&pos_reads);
        let reads = pos_reads.to_reads(true);
        for read in reads.iter() {
            println!("{}", read.to_str());
        }

        let k: usize = 8;
        let dbg_raw: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k, &reads);

        // (4) compare with true dbg
        let dbg_true_init: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k, genome);

        // (5) true k=50 (read length)
        let dbg_true: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(100, genome);

        (reads, phmm_params, dbg_raw, dbg_true_init, dbg_true)
    }

    fn generate_dataset_haploid(
        unit_seed: u64,
        hap_seed: u64,
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
        let (genome, genome_size) = generate_genome_haploid(unit_seed, hap_seed);
        println!("genome hap_a: {}", sequence_to_string(&genome[0]));

        // (2) reads and dbgs
        let (reads, phmm_params, dbg_raw, dbg_true_init, dbg_true) =
            generate_full_length_reads_and_dbgs(&genome, genome_size, read_seed);

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
        let (genome, genome_size) = generate_genome_diploid(unit_seed, hap_seed, div_seed);
        println!("genome hap_a: {}", sequence_to_string(&genome[0]));
        println!("genome hap_b: {}", sequence_to_string(&genome[1]));

        // (2) reads and dbgs
        let (reads, phmm_params, dbg_raw, dbg_true_init, dbg_true) =
            generate_full_length_reads_and_dbgs(&genome, genome_size, read_seed);

        (genome, reads, phmm_params, dbg_raw, dbg_true_init, dbg_true)
    }

    #[ignore]
    #[test]
    fn e2e_tandem_repeat_haploid() {
        let (genome, reads, phmm_params, dbg_raw, dbg_true_init, dbg_true) =
            generate_dataset_haploid(10, 0, 111);

        let scheduler = SchedulerType1::new(dbg_raw.k(), dbg_true.k(), 10.0);
        let dbg_infer = infer(&dbg_raw, &reads, &phmm_params, &scheduler, 5);

        println!("{} {}", dbg_true_init.n_traverse_choices(), dbg_true_init);
        println!("{} {}", dbg_true.n_traverse_choices(), dbg_true);
        println!("{} {}", dbg_infer.n_traverse_choices(), dbg_infer);

        let r = dbg_infer.compare(&dbg_true);
        println!("{:?}", r);
        assert_eq!(r.n_true, 408);
        assert_eq!(r.n_error, 89);

        dbg_infer.compare_with_seq(&dbg_true, &genome[0]);
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
