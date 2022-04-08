//!
//! End-to-end genome inference tests
//!
//! 1. generate genome
//! 2. generate reads from linear graph
//!
use crate::common::{sequence_to_string, Reads, Sequence};
use crate::dbg::{Dbg, HashDbg, SimpleDbg};
use crate::em::compression::compression;
use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
use crate::graph::seq_graph::SeqGraph;
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::sample::{ReadAmount, SampleProfile, StartPoints};
use crate::kmer::VecKmer;
use crate::random_seq::generate;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn e2e_compression() {
        // (1) generate genome and reads
        println!("generating genome");
        let genome = &[generate(100, 0)];
        println!("genome: {}", sequence_to_string(&genome[0]));

        println!("generating reads");
        let g = GenomeGraph::from_seqs(genome);
        let profile = ReadProfile {
            sample_profile: SampleProfile {
                read_amount: ReadAmount::Count(10),
                seed: 0,
                length: 1000,
                start_points: StartPoints::AllStartPoints,
            },
            phmm_params: PHMMParams::default(),
        };
        let reads = g.sample_reads(&profile);
        println!("reads: {}", reads.reads.len());

        // (2) crate dbg from the reads.
        //
        println!("constructing dbg");
        let k: usize = 8;
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_reads(k, &reads);
        println!("{}", dbg);

        // (3) do the inference.
        let dbg2 = compression(&dbg, &reads, &PHMMParams::default(), 1.0);
        println!("{}", dbg2);
    }
}
