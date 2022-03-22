//!
//! End-to-end genome inference tests
//!
//! 1. generate genome
//! 2. generate reads from linear graph
//!
use crate::common::{sequence_to_string, Sequence};
use crate::dbg::{Dbg, HashDbg, SimpleDbg};
use crate::em::compression::compression;
use crate::graph::genome_graph::GenomeGraph;
use crate::graph::seq_graph::SeqGraph;
use crate::hmmv2::freq::Reads;
use crate::hmmv2::params::PHMMParams;
use crate::kmer::VecKmer;
use crate::random_seq::generate;

fn sample_reads(seqs: &[Sequence], params: PHMMParams) -> Reads {
    // convert genome graph
    let g = GenomeGraph::from_seqs(seqs);
    // to (linear) phmm
    let phmm = g.to_seq_graph().to_phmm(params);
    // sample reads from the phmm
    phmm.sample_reads(100, 0, 100)
}

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
        let reads = sample_reads(genome, PHMMParams::default());
        println!("reads: {}", reads.reads.len());

        // (2) crate dbg from the reads.
        let k: usize = 8;
        let hd: HashDbg<VecKmer> = HashDbg::from_reads(k, &reads);
        println!("hashdbg {}", hd);

        let dbg = SimpleDbg::from_hashdbg(&hd);
        println!("{}", dbg);

        // (3) do the inference.
        let dbg2 = compression(&dbg, &reads, &PHMMParams::default(), depth);
    }
}
