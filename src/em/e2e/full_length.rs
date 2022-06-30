#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{sequence_to_string, Genome, Reads, Sequence};
    use crate::dbg::{Dbg, HashDbg, SimpleDbg};
    use crate::e2e::{generate_dataset, ReadType};
    use crate::em::compression::{compression, compression_step};
    use crate::em::e2e::runner::benchmark;
    use crate::em::extension::{extension, extension_step};
    use crate::em::infer;
    use crate::em::scheduler::SchedulerType1;
    use crate::genome::{simple, simple_diploid};
    use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
    use crate::graph::seq_graph::SeqGraph;
    use crate::hmmv2::params::PHMMParams;
    use crate::hmmv2::sample::{ReadAmount, SampleProfile, StartPoints};
    use crate::kmer::VecKmer;
    use crate::random_seq::generate;

    fn e2e_mock() -> (
        Genome,
        Reads,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
    ) {
        // (1) generate genome and reads
        println!("generating genome");
        let (genome, genome_size) = simple(100, 0);
        println!("genome: {}", sequence_to_string(&genome[0]));

        // (2) generate reads
        let (reads, dbg_raw, dbg_true_init, dbg_true) =
            e2e_mock_from_genome(&genome, genome_size, 10, 40);

        (genome, reads, dbg_raw, dbg_true_init, dbg_true)
    }

    fn e2e_mock_diploid() -> (
        Genome,
        Reads,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
    ) {
        // (1) generate genome and reads
        println!("generating genome");
        let (genome, genome_size) = simple_diploid();
        println!("hap1: {}", sequence_to_string(&genome[0]));
        println!("hap2: {}", sequence_to_string(&genome[1]));

        // (2) generate reads
        let (reads, dbg_raw, dbg_true_init, dbg_true) =
            e2e_mock_from_genome(&genome, genome_size, 20, 60);

        (genome, reads, dbg_raw, dbg_true_init, dbg_true)
    }

    fn e2e_mock_from_genome(
        genome: &Genome,
        genome_size: usize,
        count: usize,
        k_target: usize,
    ) -> (
        Reads,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
    ) {
        println!("generating reads");
        let dataset = generate_dataset(
            genome,
            genome_size,
            0,
            PHMMParams::default(),
            count,
            1000,
            ReadType::FullLength,
            8,
            k_target,
        );
        (
            dataset.reads,
            dataset.dbg_raw,
            dataset.dbg_true_init,
            dataset.dbg_true,
        )
    }

    #[test]
    fn e2e_compression() {
        let (genome, reads, dbg_raw, dbg_true_init, dbg_true) = e2e_mock();

        let (dbg, logs) = compression(&dbg_raw, &reads, &PHMMParams::default(), 10.0, 5);
        println!("{}", dbg);
        for log in logs.iter() {
            println!("{}", log);
        }

        let r = dbg.compare(&dbg_true_init);
        println!("{:?}", r);
        assert_eq!(r.n_true, 107);
        assert_eq!(r.n_error, 0);
    }

    #[test]
    fn e2e_extension() {
        let (genome, reads, dbg_raw, dbg_true_init, dbg_true) = e2e_mock();

        let (dbg, _) = extension(&dbg_true_init, &reads, &PHMMParams::default(), 5);
        println!("{}", dbg);
        println!("{}", dbg_true_init);
        println!("{}", dbg_true_init.n_ambiguous_intersections());
        assert_eq!(dbg.to_string(), "9,L:CCAATTCACAAAAACCACACCTTGGCCAAGGTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCACCATGGCTAGGGTCAAATCT");
    }

    #[test]
    fn e2e_full_haploid() {
        let (genome, reads, dbg_raw, dbg_true_init, dbg_true) = e2e_mock();
        println!("{}", dbg_raw.n_ambiguous_intersections());

        let (dbg_infer, r, _) = benchmark(
            &dbg_raw,
            &dbg_true,
            &reads,
            &genome,
            &PHMMParams::default(),
            10.0,
        );

        println!("{}", dbg_infer);
        println!("{}", dbg_true);
        assert_eq!(dbg_infer.to_string(), "40,L:CCAATTCACAAAAACCACACCTTGGCCAAGGTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCACCATGGCTAGGGTCAAATCT");
    }

    #[test]
    fn e2e_full_diploid() {
        let (genome, reads, dbg_raw, dbg_true_init, dbg_true) = e2e_mock_diploid();

        for read in reads.iter() {
            println!("read={}", sequence_to_string(read));
        }

        let (dbg_infer, r, _) = benchmark(
            &dbg_raw,
            &dbg_true,
            &reads,
            &genome,
            &PHMMParams::default(),
            10.0,
        );

        println!("dbg_true={}", dbg_true);
        println!("dbg_infer={}", dbg_infer);
        assert_eq!(dbg_infer.to_string(), "60,L:CCAATTCACAAAAACCACACCTTGGCCAAGCTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCATCATGGCTAGGGTCAAATCT,L:CCAATTCACAAAAACCACACCTTGGCCAAGGTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCACCATGGCTAGGGTCAAATCT");
    }
}
