#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{sequence_to_string, Genome, Reads, Sequence};
    use crate::dbg::{Dbg, HashDbg, SimpleDbg};
    use crate::e2e::{generate_dataset, Dataset, ReadType};
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

    fn e2e_mock() -> Dataset {
        // (1) generate genome and reads
        println!("generating genome");
        let (genome, genome_size) = simple(100, 0);
        println!("genome: {}", sequence_to_string(&genome[0]));

        // (2) generate reads
        e2e_mock_from_genome(&genome, genome_size, 10, 40)
    }

    fn e2e_mock_diploid() -> Dataset {
        // (1) generate genome and reads
        println!("generating genome");
        let (genome, genome_size) = simple_diploid();
        println!("hap1: {}", sequence_to_string(&genome[0]));
        println!("hap2: {}", sequence_to_string(&genome[1]));

        // (2) generate reads
        e2e_mock_from_genome(&genome, genome_size, 20, 60)
    }

    fn e2e_mock_from_genome(
        genome: &Genome,
        genome_size: usize,
        count: usize,
        k_target: usize,
    ) -> Dataset {
        println!("generating reads");
        generate_dataset(
            genome.clone(),
            genome_size,
            0,
            PHMMParams::default(),
            count,
            1000,
            ReadType::FullLength,
            8,
            k_target,
        )
    }

    #[test]
    fn e2e_compression() {
        let dataset = e2e_mock();

        let (dbg, logs) = compression(
            &dataset.dbg_raw,
            &dataset.reads,
            &PHMMParams::default(),
            10.0,
            5,
        );
        println!("{}", dbg);
        for log in logs.iter() {
            println!("{}", log);
        }

        let r = dbg.compare(&dataset.dbg_true_init);
        println!("{:?}", r);
        assert_eq!(r.n_true, 107);
        assert_eq!(r.n_error, 0);
    }

    #[test]
    fn e2e_extension() {
        let dataset = e2e_mock();

        let (dbg, _) = extension(
            &dataset.dbg_true_init,
            &dataset.reads,
            &PHMMParams::default(),
            5,
        );
        println!("{}", dbg);
        println!("{}", dataset.dbg_true_init);
        println!("{}", dataset.dbg_true_init.n_ambiguous_intersections());
        assert_eq!(dbg.to_string(), "9,L:CCAATTCACAAAAACCACACCTTGGCCAAGGTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCACCATGGCTAGGGTCAAATCT");
    }

    #[test]
    fn e2e_full_haploid() {
        let dataset = e2e_mock();
        println!("{}", dataset.dbg_raw.n_ambiguous_intersections());

        let (dbg_infer, r, _) = benchmark(&dataset, 10.0);

        println!("{}", dbg_infer);
        println!("{}", dataset.dbg_true);
        assert_eq!(dbg_infer.to_string(), "40,L:CCAATTCACAAAAACCACACCTTGGCCAAGGTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCACCATGGCTAGGGTCAAATCT");
    }

    #[test]
    fn e2e_full_diploid() {
        let dataset = e2e_mock_diploid();

        for read in dataset.reads.iter() {
            println!("read={}", sequence_to_string(read));
        }

        let (dbg_infer, r, _) = benchmark(&dataset, 10.0);

        println!("dbg_true={}", dataset.dbg_true);
        println!("dbg_infer={}", dbg_infer);
        assert_eq!(dbg_infer.to_string(), "60,L:CCAATTCACAAAAACCACACCTTGGCCAAGCTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCATCATGGCTAGGGTCAAATCT,L:CCAATTCACAAAAACCACACCTTGGCCAAGGTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCACCATGGCTAGGGTCAAATCT");
    }
}
