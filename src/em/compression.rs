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
use crate::min_flow::min_cost_flow_convex_fast;

///
/// Compression algorithm
///
/// ## TODOs
///
/// * avoid dbg copy?
///
pub fn compression<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    depth: Freq,
) -> Dbg<N, E> {
    // e-step
    // calculate node_freqs by using current dbg.
    let node_freqs = compression_e_step(dbg, reads, params);
    println!("node_freqs={}", node_freqs);

    // m-step
    // convert it to the
    let copy_nums = compression_m_step(dbg, &node_freqs, depth);
    println!("copy_nums={}", copy_nums);

    let mut new_dbg = dbg.clone();
    new_dbg.set_node_copy_nums(&copy_nums);
    new_dbg
}

///
/// E-step of compression
///
/// ## Details
///
/// * convert dbg into phmm.
/// * calculate node frequencies by forward/backward algorithm to emit the reads.
///
fn compression_e_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
) -> NodeFreqs {
    let phmm = dbg.to_phmm(params.clone());
    phmm.to_node_freqs(reads)
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
fn compression_m_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    node_freqs: &NodeFreqs,
    depth: Freq,
) -> NodeCopyNums {
    let node_freqs = node_freqs.clone() / depth;
    let edbg = dbg.to_edbg_with_attr(Some(&node_freqs));
    let flow = min_cost_flow_convex_fast(&edbg.graph);
    match flow {
        None => panic!(),
        // an edge in edbg corresponds to a node in dbg
        // so edgevec for edbg can be converted to nodevec for dbg.
        Some(copy_nums) => copy_nums.switch_index(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks::*;

    #[test]
    fn compression_m_step_dbg() {
        let dbg = mock_base();
        let node_freqs = NodeFreqs::new(dbg.n_nodes(), 1.9);
        let copy_nums = compression_m_step(&dbg, &node_freqs, 1.0);
        println!("{}", copy_nums);
        assert_eq!(copy_nums.to_vec(), vec![2; dbg.n_nodes()]);

        let copy_nums = compression_m_step(&dbg, &node_freqs, 2.0);
        println!("{}", copy_nums);
        assert_eq!(copy_nums.to_vec(), vec![1; dbg.n_nodes()]);
    }

    #[test]
    fn compression_dbg() {
        let dbg = mock_intersection();
        let reads = Reads {
            reads: vec![
                b"AACTAGCTT".to_vec(),
                b"AACTAGCTT".to_vec(),
                b"AACTAGCTT".to_vec(),
            ],
        };
        let params = PHMMParams::default();
        let dbg = compression(&dbg, &reads, &params, 3.0);
    }
}
