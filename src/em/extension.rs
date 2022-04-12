//!
//! **Extension**
//!
//! convert k-dbg into k+1-dbg.
//!
//! ## E-step
//!
//!
//! ## M-step
//!
//!
use crate::common::{CopyNum, Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, EdgeCopyNums};
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::trans_table::EdgeFreqs;
use crate::kmer::kmer::KmerLike;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
pub mod flow_intersection;
use flow_intersection::{FlowIntersection, FlowIntersectionEdge, FlowIntersectionNode};
pub mod intersection_graph;

///
/// Extension algorithm
///
/// ## Details
///
/// * e-step
///     Calculate edge_freqs on dbg.
///
/// * m-step
///     Maximize the score for each intersections.
///
/// ## TODOs
///
/// * avoid dbg copy?
///
pub fn extension<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
) -> Dbg<N, E> {
    // (1) e-step infer edge freqs
    println!("extension::e_step");
    let edge_freqs = e_step(dbg, reads, params);
    println!("edge_freqs={}", edge_freqs);

    // (2) m-step infer the best copy nums
    println!("extension::m_step");
    let copy_nums = m_step(dbg, &edge_freqs);
    println!("copy_nums={}", copy_nums);

    let mut new_dbg = dbg.clone();
    new_dbg.set_edge_copy_nums(Some(&copy_nums));
    new_dbg
}

///
/// E-step of extension
///
/// ## Details
///
/// * convert dbg into phmm.
/// * calculate edge frequencies by forward/backward algorithm to emit the reads.
///
fn e_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
) -> EdgeFreqs {
    let phmm = dbg.to_phmm(params.clone());
    phmm.to_edge_freqs(reads)
}

///
/// M-step of extension
///
/// ## Params
///
/// * `dbg`
///     de bruijn graph
/// * `edge_freqs`
///     edge usage frequencies.
///
/// ## Details
///
///
fn m_step<N: DbgNode, E: DbgEdge>(dbg: &Dbg<N, E>, edge_freqs: &EdgeFreqs) -> EdgeCopyNums {
    let default_value = 0;
    let mut ecn = EdgeCopyNums::new(dbg.n_edges(), default_value);
    for fi in dbg.iter_flow_intersections(edge_freqs) {
        // get an optimized flow intersection
        let fio = fi.convert();

        println!("extension iter m {} {}", fi, fio);

        // check if there is no inconsistent edge copy numbers.
        assert!(fio.has_valid_node_copy_nums());
        assert!(fio.all_edges_has_copy_num());

        // store fio's edge copy number information into ecn vector.
        for (_, _, e) in fio.bi.iter_edges() {
            assert!(ecn[e.index] == default_value);
            ecn[e.index] = e.copy_num.unwrap();
        }
    }
    ecn
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks::{mock_base, mock_intersection};
    use crate::hmmv2::trans_table::EdgeFreqs;

    #[test]
    fn em_extension_m_mock_base() {
        let dbg = mock_base();
        let freqs = EdgeFreqs::new(dbg.n_edges(), 1.1);
        println!("{}", dbg);
        println!("{}", freqs);
        let copy_nums = m_step(&dbg, &freqs);
        println!("{}", copy_nums);
        assert_eq!(copy_nums.to_vec(), vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);
    }

    #[test]
    fn em_extension_e_mock_intersection() {
        let dbg = mock_intersection();
        let read = b"AACTAGCTT";
        let reads = Reads {
            reads: vec![read.to_vec()],
        };
        println!("{}", dbg);
        let params = PHMMParams::default();

        let freqs = e_step(&dbg, &reads, &params);
        let freqs_true = EdgeFreqs::from_slice(
            &[
                0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0,
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
            ],
            0.0,
        );
        println!("{}", freqs);
        println!("{}", freqs.diff_f64(&freqs_true));
        assert!(freqs.diff_f64(&freqs_true) < 0.05);

        // TODO allow fragmented seq in to_copy_nums_of_seq
        let (_, ecn_true) = dbg.to_copy_nums_of_seq(read).unwrap();
        println!("{}", ecn_true);
    }

    #[test]
    fn em_extension_m_mock_intersection() {
        let dbg = mock_intersection();
        println!("{}", dbg);
        let freqs = EdgeFreqs::from_slice(
            &[
                0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0,
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
            ],
            0.0,
        );
        let copy_nums = m_step(&dbg, &freqs);
        println!("{}", copy_nums);
        // assert_eq!(copy_nums.to_vec(), vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);
    }
}
