#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{sequence_to_string, Genome, Reads, Sequence};
    use crate::dbg::{Dbg, HashDbg, SimpleDbg};
    use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
    use crate::graph::seq_graph::SeqGraph;
    use crate::hmmv2::params::PHMMParams;
    use crate::hmmv2::sample::{ReadAmount, SampleProfile, StartPoints};
    use crate::kmer::VecKmer;
    use crate::random_seq::generate;

    fn generate_e2e_fragment_mock() {
        // (1) generate genome and reads
        println!("generating genome");
        let genome = vec![generate(500, 0)];
        println!("genome: {}", sequence_to_string(&genome[0]));

        println!("generating reads");
        let g = GenomeGraph::from_seqs(&genome);
        let profile = ReadProfile {
            sample_profile: SampleProfile {
                read_amount: ReadAmount::TotalBases(100),
                seed: 0,
                length: 100,
                start_points: StartPoints::Random,
                endable: false,
            },
            phmm_params: PHMMParams::default(),
        };
        let reads = g.sample_reads(&profile);
        println!("n_reads: {}", reads.len());
        assert_eq!(reads.len(), 2);
        for (i, read) in reads.iter().enumerate() {
            println!("{}", sequence_to_string(read));
            if i == 0 {
                assert_eq!(sequence_to_string(read), "TGAATCCTAGATCCCGTTGTCGGGGCTCGGCGTTTGCTTTCTTAGATTCCGATAAGTAGATGGTTTCCTGGGTGAGGGCACTATTAAAGCGGCGATTTG");
            } else if i == 1 {
                assert_eq!(sequence_to_string(read), "AGCGATTAAACACCCTATAAAAATGGCCATCCGCTGAGCTTGCATCACAGTTGGTCTTACACATGCCTGCTTCATCAAAGTCCCACTGCGCCATCA");
            }
        }
    }

    #[test]
    fn e2e_fragment() {
        generate_e2e_fragment_mock();
    }
}
