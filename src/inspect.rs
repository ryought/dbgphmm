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
use crate::e2e::{
    generate_500bp_case_dataset, generate_difficult_diploid_tandem_repeat_dataset_full_length,
    Dataset,
};
use crate::hmmv2::common::{PHMMEdge, PHMMNode};
use crate::kmer::veckmer::kmer;
use crate::kmer::VecKmer;
use petgraph::graph::{EdgeIndex, NodeIndex};

///
///
///
pub fn generate_500bp_case() -> (Dataset, SimpleDbg<VecKmer>, SimpleDbg<VecKmer>) {
    let dataset = generate_500bp_case_dataset();
    let k = 12;

    // true dbg
    let dbg_true = SimpleDbg::from_styled_seqs(k, dataset.genome());
    assert!(dbg_true.is_valid());

    // false MAP dbg
    let mut dbg_opt = dbg_true.clone();
    dbg_opt.edit_copy_nums_by_seq(b"GACAAGCTAGGCAAGCTAGGACAAGCTAGG", -1);
    // dbg_opt.edit_copy_nums_by_seq(b"GACAAGCTAGGACAAGCTAGG", -1);
    assert!(dbg_opt.is_valid());

    (dataset, dbg_true, dbg_opt)
}

///
///
///
pub fn generate_case3() -> (
    Dataset,
    SimpleDbg<VecKmer>,
    SimpleDbg<VecKmer>,
    Vec<EdgeIndex>,
) {
    let dataset = generate_difficult_diploid_tandem_repeat_dataset_full_length();
    let k = 12;

    // true dbg
    let dbg_true = SimpleDbg::from_styled_seqs(k, dataset.genome());

    // false MAP dbg
    let mut dbg_opt = dbg_true.clone();
    dbg_opt.edit_copy_nums_by_seq(b"AAACCACACCTGGCCAAGGTATC", 1);
    dbg_opt.edit_copy_nums_by_seq(b"AAACCACACCTTTGCCAAGGTATC", -1);

    // edges
    let mut edges = Vec::new();
    // check trans_prob difference
    let phmm_true = dbg_true.to_phmm(dataset.params());
    let phmm_opt = dbg_opt.to_phmm(dataset.params());
    // init_prob of almost all nodes will be changed
    // because the total genome_size was changed.
    // for (node, _) in dbg_true.nodes() {
    //     let pi_t = phmm_true.node(node).init_prob();
    //     let pi_o = phmm_opt.node(node).init_prob();
    //     if pi_t != pi_o {
    //         println!("n{} {} {}", node.index(), pi_t, pi_o);
    //     }
    // }
    for (edge, s, t, _) in dbg_true.edges() {
        let pt_t = phmm_true.edge(edge).trans_prob();
        let pt_o = phmm_opt.edge(edge).trans_prob();
        if pt_t != pt_o {
            println!(
                "e{} {}(x{},x{}) -> {}(x{},x{}) {} {}",
                edge.index(),
                dbg_true.kmer(s),
                dbg_true.copy_num(s),
                dbg_opt.copy_num(s),
                dbg_true.kmer(t),
                dbg_true.copy_num(t),
                dbg_opt.copy_num(t),
                pt_t,
                pt_o
            );
            edges.push(edge);
        }
    }

    (dataset, dbg_true, dbg_opt, edges)
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn e2e_difficult_tandem_repeat_generation() {
        let (dataset, dbg_true, dbg_opt, _) = generate_case3();
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

    #[ignore]
    #[test]
    fn e2e_difficult_tandem_repeat_show_reads() {
        let (dataset, dbg_true, dbg_opt, _) = generate_case3();
        dataset.show_reads_with_genome();
    }

    #[ignore]
    #[test]
    fn e2e_difficult_tandem_repeat_compare_mapping() {
        let (dataset, dbg_true, dbg_opt, _) = generate_case3();
        let copy_nums_true = dbg_true.to_node_copy_nums();
        let copy_nums_opt = dbg_opt.to_node_copy_nums();

        dbg_true.compare_mappings(
            dataset.params(),
            dataset.reads(),
            &copy_nums_true,
            &copy_nums_opt,
        );
    }

    #[test]
    fn e2e_500bp_generation() {
        let (dataset, dbg_true, dbg_opt) = generate_500bp_case();
        dataset.show_reads_with_genome();
        println!("{}", dbg_true);
        println!("{}", dbg_opt);
        let p_true = dbg_true.to_full_prob(dataset.params(), dataset.reads());
        println!("p_true={}", p_true);
        let p_opt = dbg_opt.to_full_prob(dataset.params(), dataset.reads());
        println!("p_opt={}", p_opt);
        let missing_node = dbg_opt.find_node_from_kmer(&kmer(b"AAGCTAGGCAAG")).unwrap();
        assert_eq!(dbg_true.copy_num(missing_node), 1);
        assert_eq!(dbg_opt.copy_num(missing_node), 0);
    }
}
