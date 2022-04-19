//!
//! Compression
//!
//! ## E-step
//!
//! Estimate node freq.
//!
//! ## M-step
//!
//! Solve min-flow and determine node copy number.
//!
use crate::common::{Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::hmmv2::freq::NodeFreqs;
use crate::hmmv2::params::PHMMParams;
use crate::min_flow::{min_cost_flow_convex_fast, total_cost, Cost};
use crate::prob::Prob;

///
/// Log information store of each iteration in compression
///
#[derive(Clone, Debug)]
pub struct CompressionLog {
    /// Full probability
    full_prob: Prob,
    /// Min-flow error
    min_flow_score: Cost,
}

impl CompressionLog {
    pub fn new(full_prob: Prob, min_flow_score: Cost) -> Self {
        CompressionLog {
            full_prob,
            min_flow_score,
        }
    }
}

///
/// Compression full algorithm by running `compression_step` iteratively.
///
/// * max_iter: max iteration loop count of EM.
///
pub fn compression<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    depth: Freq,
    max_iter: usize,
) -> (Dbg<N, E>, Vec<CompressionLog>) {
    let depths = vec![depth; max_iter];
    compression_with_depths(dbg, reads, params, &depths)
}

///
/// Compression full algorithm by running `compression_step` iteratively
/// with specified depths.
///
pub fn compression_with_depths<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    depths: &[Freq],
) -> (Dbg<N, E>, Vec<CompressionLog>) {
    let mut dbg = dbg.clone();
    let mut logs = Vec::new();

    // iterate EM steps
    for &depth in depths {
        let (dbg_new, is_updated, log) = compression_step(&dbg, reads, params, depth);
        logs.push(log);

        // if the single EM step does not change the DBG model, stop iteration.
        if !is_updated {
            break;
        }
        dbg = dbg_new;
    }

    (dbg, logs)
}

///
/// Compression algorithm
///
/// ## TODOs
///
/// * avoid dbg copy?
///
pub fn compression_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    depth: Freq,
) -> (Dbg<N, E>, bool, CompressionLog) {
    // e-step
    // calculate node_freqs by using current dbg.
    let (node_freqs, full_prob) = e_step(dbg, reads, params);

    // m-step
    let (copy_nums, min_flow_score) = m_step(dbg, &node_freqs, depth);

    let mut new_dbg = dbg.clone();
    let is_updated = new_dbg.set_node_copy_nums(&copy_nums);
    new_dbg.remove_zero_copy_node();

    // create log
    let log = CompressionLog::new(full_prob, min_flow_score);
    println!("{}", new_dbg);

    (new_dbg, is_updated, log)
}

///
/// E-step of compression
///
/// ## Details
///
/// * convert dbg into phmm.
/// * calculate node frequencies by forward/backward algorithm to emit the reads.
///
fn e_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
) -> (NodeFreqs, Prob) {
    let phmm = dbg.to_phmm(params.clone());
    phmm.to_node_freqs_parallel(reads)
}

///
/// M-step of compression
///
/// ## Params
///
/// * `dbg`
///     de bruijn graph
/// * `node_freqs`
///     node usage frequencies.
/// * `depth`
///     expected depth of each node.
///
/// ## Details
///
/// * convert into edge-centric dbg (with freq)
/// * solve min-flow with demand=0 capacity=infty cost=squared-error-from-freq
///
fn m_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    node_freqs: &NodeFreqs,
    depth: Freq,
) -> (NodeCopyNums, Cost) {
    let node_freqs = node_freqs.clone() / depth;
    let edbg = dbg.to_edbg_with_attr(Some(&node_freqs));
    let flow = min_cost_flow_convex_fast(&edbg.graph);
    match flow {
        None => panic!("compression::m_step cannot find optimal flow."),
        // an edge in edbg corresponds to a node in dbg
        // so edgevec for edbg can be converted to nodevec for dbg.
        Some(copy_nums) => {
            let cost = total_cost(&edbg.graph, &copy_nums);
            (copy_nums.switch_index(), cost)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks::*;

    #[test]
    fn em_compression_m_step_dbg() {
        let dbg = mock_base();
        let node_freqs = NodeFreqs::new(dbg.n_nodes(), 1.9);
        let (copy_nums, score) = m_step(&dbg, &node_freqs, 1.0);
        println!("{}", copy_nums);
        println!("score={}", score);
        assert_eq!(copy_nums.to_vec(), vec![2; dbg.n_nodes()]);
        assert_abs_diff_eq!(score, 0.07);

        let (copy_nums, score) = m_step(&dbg, &node_freqs, 2.0);
        println!("{}", copy_nums);
        println!("{}", node_freqs);
        println!("score={}", score);
        assert_abs_diff_eq!(score, 0.0175);
        assert_eq!(copy_nums.to_vec(), vec![1; dbg.n_nodes()]);
    }

    #[test]
    fn em_compression_step_intersection() {
        let dbg = mock_intersection();
        let reads = Reads {
            reads: vec![
                b"AACTAGCTT".to_vec(),
                b"AACTAGCTT".to_vec(),
                b"AACTAGCTT".to_vec(),
            ],
        };
        let params = PHMMParams::default();
        println!("{}", dbg);
        assert_eq!(dbg.genome_size(), 18);
        assert_eq!(dbg.to_string(), "4,L:AACTAGCTT,L:CCGTAGGGC");

        let (dbg_v2, is_updated, _) = compression_step(&dbg, &reads, &params, 3.0);
        println!("{}", dbg_v2);
        println!("{}", dbg_v2.genome_size());
        println!("is_updated={}", is_updated);
        assert_eq!(dbg_v2.genome_size(), 9);
        assert_eq!(dbg_v2.to_string(), "4,L:AACTAGCTT");
        assert!(is_updated);

        // compress again
        let (dbg_v3, is_updated, _) = compression_step(&dbg_v2, &reads, &params, 3.0);
        println!("{}", dbg_v3);
        println!("{}", dbg_v3.genome_size());
        println!("is_updated={}", is_updated);
        assert_eq!(dbg_v3.genome_size(), 9);
        assert_eq!(dbg_v3.to_string(), "4,L:AACTAGCTT");
        assert!(!is_updated);
    }
}
