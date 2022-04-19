#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{sequence_to_string, Genome, Reads, Sequence};
    use crate::dbg::{Dbg, HashDbg, SimpleDbg};
    use crate::em::compression::{compression, compression_step, compression_with_depths};
    use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
    use crate::graph::seq_graph::SeqGraph;
    use crate::hmmv2::params::PHMMParams;
    use crate::hmmv2::sample::{ReadAmount, SampleProfile, StartPoints};
    use crate::kmer::VecKmer;
    use crate::random_seq::generate;

    fn generate_e2e_fragment_mock() -> (Genome, Reads, SimpleDbg<VecKmer>, SimpleDbg<VecKmer>) {
        // (1) generate genome and reads
        println!("generating genome");
        let genome_size = 200;
        let genome = vec![generate(genome_size, 0)];
        println!("genome: {}", sequence_to_string(&genome[0]));

        println!("generating reads");
        let g = GenomeGraph::from_seqs(&genome);
        let profile = ReadProfile {
            sample_profile: SampleProfile {
                read_amount: ReadAmount::TotalBases(genome_size * 10),
                seed: 0,
                length: 50,
                start_points: StartPoints::Random,
                endable: false,
            },
            phmm_params: PHMMParams::default(),
        };
        let (reads, pos) = g.sample_reads_with_pos(&profile);
        for (i, read) in reads.iter().enumerate() {
            println!("{}", sequence_to_string(read));
            println!("{:?}", pos[i]);
        }

        let k: usize = 8;
        let dbg_raw: SimpleDbg<VecKmer> = SimpleDbg::from_reads(k, &reads);
        println!("{}", dbg_raw);

        // (4) compare with true dbg
        let dbg_true: SimpleDbg<VecKmer> = SimpleDbg::from_seq(k, &genome[0]);

        (genome, reads, dbg_raw, dbg_true)
    }

    #[test]
    fn e2e_fragment() {
        let (genome, reads, dbg_raw, dbg_true) = generate_e2e_fragment_mock();

        // let (dbg, _) = compression(&dbg_raw, &reads, &PHMMParams::default(), 1.0, 10);
        let (dbg, logs) = compression_with_depths(
            &dbg_raw,
            &reads,
            &PHMMParams::default(),
            &[
                1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 5.0, 5.0, 8.0, 8.0, 10.0, 10.0,
            ],
        );
        println!("{}", dbg);
        println!("{:?}", logs);
        println!("{}", dbg_true);

        // let r = dbg.compare(&dbg_true);
        // println!("{:?}", r);
    }
}
