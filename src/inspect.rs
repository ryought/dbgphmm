//!
//! inspection of inference failure cases
//!
//! * tandem repeat
//!
use crate::common::{
    sequence_to_string, Genome, PositionedReads, PositionedSequence, Reads, Seq, Sequence,
    StyledSequence,
};
use crate::dbg::{Dbg, HashDbg, SimpleDbg};
use crate::e2e::{generate_difficult_diploid_tandem_repeat_dataset_full_length, Dataset};
use crate::kmer::veckmer::kmer;
use crate::kmer::VecKmer;

///
///
///
pub fn generate_case3() -> (Dataset, SimpleDbg<VecKmer>, SimpleDbg<VecKmer>) {
    let dataset = generate_difficult_diploid_tandem_repeat_dataset_full_length();
    let k = 12;

    // true dbg
    let dbg_true = SimpleDbg::from_styled_seqs(k, dataset.genome());

    // false MAP dbg
    let mut copy_nums_opt = dbg_true.to_node_copy_nums();
    for node in dbg_true
        .to_nodes_of_styled_seq(&StyledSequence::linear_fragment(
            b"AAACCACACCTGGCCAAGGTATC".to_vec(),
        ))
        .unwrap()
    {
        copy_nums_opt[node] += 1;
    }
    for node in dbg_true
        .to_nodes_of_styled_seq(&StyledSequence::linear_fragment(
            b"AAACCACACCTTTGCCAAGGTATC".to_vec(),
        ))
        .unwrap()
    {
        copy_nums_opt[node] -= 1;
    }
    let mut dbg_opt = dbg_true.clone();
    dbg_opt.set_node_copy_nums(&copy_nums_opt);

    (dataset, dbg_true, dbg_opt)
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn e2e_difficult_tandem_repeat_generation() {
        let (dataset, dbg_true, dbg_opt) = generate_case3();
        assert_eq!(dbg_true.genome_size(), 2209);
        assert_eq!(dbg_opt.genome_size(), 2208);

        let missing_node = dbg_opt.find_node_from_kmer(&kmer(b"ACACCTTTGCCA")).unwrap();
        assert_eq!(dbg_true.copy_num(missing_node), 1);
        assert_eq!(dbg_opt.copy_num(missing_node), 0);

        let p_true = dbg_true.to_full_prob(dataset.params(), dataset.reads());
        println!("p_true={}", p_true);
        assert!((p_true.to_log_value() - (-17125.0)).abs() < 1.0);
        let p_opt = dbg_opt.to_full_prob(dataset.params(), dataset.reads());
        println!("p_opt={}", p_opt);
        assert!((p_opt.to_log_value() - (-17121.0)).abs() < 1.0);
    }
}
