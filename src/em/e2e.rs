//!
//! End-to-end genome inference tests
//!
//! It starts from a raw dbg constructed from reads.
//!
//! 1. generate genome
//! 2. generate reads from linear graph
//!
use crate::common::{sequence_to_string, Genome, Reads, Sequence};
use crate::dbg::{Dbg, HashDbg, SimpleDbg};
use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
use crate::graph::seq_graph::SeqGraph;
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::sample::{ReadAmount, SampleProfile, StartPoints};
use crate::kmer::VecKmer;
use crate::random_seq::generate;
mod fragments;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::em::compression::{compression, compression_step};
    use crate::em::extension::{extension, extension_step};

    fn e2e_mock() -> (Genome, Reads, SimpleDbg<VecKmer>, SimpleDbg<VecKmer>) {
        // (1) generate genome and reads
        println!("generating genome");
        let genome = vec![generate(100, 0)];
        println!("genome: {}", sequence_to_string(&genome[0]));

        println!("generating reads");
        let g = GenomeGraph::from_seqs(&genome);
        let profile = ReadProfile {
            sample_profile: SampleProfile {
                read_amount: ReadAmount::Count(10),
                seed: 0,
                length: 1000,
                start_points: StartPoints::AllStartPoints,
                endable: false,
            },
            phmm_params: PHMMParams::default(),
        };
        let reads = g.sample_reads(&profile);
        println!("n_reads: {}", reads.len());

        // (2) crate dbg from the reads.
        //
        println!("constructing dbg");
        let k: usize = 8;
        let dbg_raw: SimpleDbg<VecKmer> = SimpleDbg::from_reads(k, &reads);
        println!("{}", dbg_raw);

        // (4) compare with true dbg
        let dbg_true: SimpleDbg<VecKmer> = SimpleDbg::from_seq(k, &genome[0]);

        (genome, reads, dbg_raw, dbg_true)
    }

    #[test]
    fn e2e_compression() {
        let (genome, reads, dbg_raw, dbg_true) = e2e_mock();

        let (dbg, logs) = compression(&dbg_raw, &reads, &PHMMParams::default(), 10.0, 5);
        println!("{}", dbg);
        println!("{:?}", logs);

        let r = dbg.compare(&dbg_true);
        println!("{:?}", r);
        assert_eq!(r.n_true, 107);
        assert_eq!(r.n_error, 0);
    }

    #[test]
    fn e2e_extension() {
        let (genome, reads, dbg_raw, dbg_true) = e2e_mock();

        let (dbg, _) = extension(&dbg_raw, &reads, &PHMMParams::default(), 5);
        println!("{}", dbg);
    }
}
