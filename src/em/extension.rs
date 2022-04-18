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
use crate::prob::Prob;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
pub mod flow_intersection;
use flow_intersection::{FlowIntersection, FlowIntersectionEdge, FlowIntersectionNode};
pub mod intersection_graph;

///
/// Log information store of each iteration in extension
///
pub struct ExtensionLog {
    /// Full probability
    full_prob: Prob,
}

impl ExtensionLog {
    pub fn new(full_prob: Prob) -> Self {
        ExtensionLog { full_prob }
    }
}

///
/// Extension full algorithm by running `extension_step()` iteratively.
///
/// convert k-dBG into k+1-dBG.
///
/// * max_iter: max iteration loop count
///
pub fn extension<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    max_iter: usize,
) -> (Dbg<N, E>, Vec<ExtensionLog>) {
    let mut dbg = dbg.clone();
    let mut logs = Vec::new();

    // iterate EM steps
    for i in 0..max_iter {
        println!("extension {}th iteration", i);
        let (dbg_new, is_updated, log) = extension_step(&dbg, reads, params);
        logs.push(log);
        if !is_updated {
            break;
        }
        dbg = dbg_new;
    }

    // convert to k+1 dbg
    (dbg.to_kp1_dbg(), logs)
}

///
/// Extension algorithm
///
/// It returns
/// * the updated dBG
/// * the dBG was changed or not.
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
pub fn extension_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
) -> (Dbg<N, E>, bool, ExtensionLog) {
    // (1) e-step infer edge freqs
    println!("extension::e_step");
    let (edge_freqs, full_prob) = e_step(dbg, reads, params);
    println!("edge_freqs={}", edge_freqs);
    println!("full_prob={}", full_prob);

    // (2) m-step infer the best copy nums
    println!("extension::m_step");
    let copy_nums = m_step(dbg, &edge_freqs);
    println!("copy_nums={}", copy_nums);

    let mut new_dbg = dbg.clone();
    let is_updated = new_dbg.set_edge_copy_nums(Some(&copy_nums));
    println!("is_updated={}", is_updated);

    (new_dbg, is_updated, ExtensionLog::new(full_prob))
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
) -> (EdgeFreqs, Prob) {
    let phmm = dbg.to_phmm(params.clone());
    phmm.to_edge_freqs_parallel(reads)
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
        if !fi.is_tip_intersection() {
            // get an optimized flow intersection
            let (fio, cost) = fi.convert();

            if !fi.can_uniquely_convertable() {
                println!("extension optimized iter m {} {}", fi, fio);
            }

            // check if there is no inconsistent edge copy numbers.
            assert!(fio.has_valid_node_copy_nums());
            assert!(fio.is_resolved());

            // store fio's edge copy number information into ecn vector.
            for (_, _, e) in fio.bi.iter_edges() {
                assert!(ecn[e.index] == default_value);
                ecn[e.index] = e.copy_num.unwrap();
            }
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
        assert_eq!(copy_nums.to_vec(), vec![1, 1, 1, 1, 1, 1, 0, 1, 1, 1]);
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

        let (freqs, _) = e_step(&dbg, &reads, &params);
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

    #[test]
    fn em_extension_step_intersection() {
        let dbg = mock_intersection();
        let read = b"AACTAGCTT";
        let reads = Reads {
            reads: vec![read.to_vec()],
        };
        println!("{}", dbg);
        let params = PHMMParams::default();

        let (dbg_v2, is_updated, _) = extension_step(&dbg, &reads, &params);
        println!("{}", dbg_v2);
        println!("is_updated={}", is_updated);
        assert!(is_updated);

        // extension again
        let (dbg_v3, is_updated, _) = extension_step(&dbg_v2, &reads, &params);
        println!("{}", dbg_v2);
        println!("is_updated={}", is_updated);
        assert!(!is_updated);
    }

    #[test]
    fn em_extension_all_intersection() {
        let dbg = mock_intersection();
        let read = b"AACTAGCTT";
        let reads = Reads {
            reads: vec![read.to_vec(), read.to_vec(), read.to_vec()],
        };
        println!("{}", dbg);
        assert_eq!(format!("{}", dbg), "4,L:AACTAGGGC,L:CCGTAGCTT");
        println!("genome_size={}", dbg.genome_size());
        assert_eq!(dbg.genome_size(), 18);
        let params = PHMMParams::default();

        // loop
        let (dbg_extended, _) = extension(&dbg, &reads, &params, 5);
        println!("{}", dbg_extended);
        assert_eq!(format!("{}", dbg_extended), "5,L:AACTAGCTT,L:CCGTAGGGC");
        println!("genome_size={}", dbg_extended.genome_size());
        assert_eq!(dbg_extended.genome_size(), 18);
    }
}
