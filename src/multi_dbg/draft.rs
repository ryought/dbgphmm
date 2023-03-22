//!
//! Constructor of draft dbg from reads or genomes
//!
use super::MultiDbg;
use crate::common::collection::starts_and_ends_of_genome;
use crate::common::{Seq, StyledSequence};
use crate::dbg::{draft::EndNodeInference, Dbg, SimpleDbg};
use crate::e2e::Dataset;
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
        end_node_inference: &EndNodeInference<VecKmer>,
    ) -> Self
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        let dbg: SimpleDbg<VecKmer> =
            SimpleDbg::create_draft_from_fragment_seqs_with_adjusted_coverage(
                k,
                seqs,
                base_coverage,
                ave_read_length,
                p_error,
                end_node_inference,
            );
        dbg.into()
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
    ///
    /// Create read-draft MultiDbg from dataset
    ///
    /// * use end inference from genome
    /// * true path
    ///
    pub fn create_draft_from_dataset(k: usize, dataset: &Dataset) -> Self {
        let d = Self::create_draft_from_reads(
            k,
            dataset.reads(),
            dataset.coverage(),
            dataset.average_read_length(),
            dataset.params().p_error().to_value(),
            // &EndNodeInference::Auto,
            &EndNodeInference::Custom(starts_and_ends_of_genome(dataset.genome(), k)),
        );
        d
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
}
