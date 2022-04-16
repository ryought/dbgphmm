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
                read_amount: ReadAmount::Count(5),
                seed: 0,
                length: 100,
                start_points: StartPoints::Random,
            },
            phmm_params: PHMMParams::default(),
        };
        let reads = g.sample_reads(&profile);
        println!("n_reads: {}", reads.len());
    }

    #[test]
    fn e2e_fragment() {
        generate_e2e_fragment_mock();
    }
}
