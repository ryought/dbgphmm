//!
//! Constructor of draft dbg from reads or genomes
//!
use super::MultiDbg;
use crate::common::{Seq, StyledSequence};
use crate::dbg::{Dbg, SimpleDbg};
use crate::kmer::VecKmer;
use crate::utils::timer;

impl MultiDbg {
    ///
    /// Create from reads
    ///
    pub fn create_draft_from_reads<T>(
        k: usize,
        seqs: T,
        base_coverage: f64,
        ave_read_length: usize,
        p_error: f64,
    ) -> Self
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        // let dbg: SimpleDbg<VecKmer> =
        //     SimpleDbg::create_draft_from_fragment_seqs_with_adjusted_coverage(
        //         k,
        //         seqs,
        //         base_coverage,
        //         ave_read_length,
        //         p_error,
        //     );
        // dbg.into()
        unimplemented!();
    }
    ///
    /// Create from styled seqs
    ///
    pub fn create_from_styled_seqs<T>(k: usize, seqs: T) -> Self
    where
        T: IntoIterator,
        T::Item: AsRef<StyledSequence>,
    {
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_styled_seqs(k, seqs);
        dbg.into()
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome;

    #[test]
    #[ignore]
    fn from_styled_seqs_large() {
        let (genome, genome_size) = genome::tandem_repeat_polyploid_with_unique_homo_ends(
            1_000, 1_000, 0, 1_000, 2, 0.01, 0,
        );
        let (mdbg, t) = timer(|| MultiDbg::create_from_styled_seqs(40, &genome));
        // ~2391ms
        println!("created mdbg in {}ms", t);
        mdbg.to_gfa_file("g1m.gfa");
        let (m, t) = timer(|| mdbg.to_kmer_map());
        // ~182ms
        println!("created map in {}ms", t);
        println!("m {}", m.len());
    }

    #[test]
    fn from_reads() {}
}
