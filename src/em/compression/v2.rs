//!
//! CompressionV2
//!
use crate::common::{CopyNum, Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::dbg::edge_centric::impls::{SimpleEDbgEdgeWithAttr, MAX_COPY_NUM_OF_EDGE};
use crate::hmmv2::freq::NodeFreqs;
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::trans_table::EdgeFreqs;
use crate::kmer::kmer::{Kmer, KmerLike};
use crate::min_flow::convex::ConvexCost;
use crate::min_flow::flow::FlowEdge;
use crate::min_flow::{min_cost_flow_convex_fast, total_cost, Cost};
use crate::prob::Prob;
use crate::vector::{DenseStorage, EdgeVec, NodeVec, Storage};

#[derive(Clone, Debug, Copy, PartialEq, Default)]
pub struct CompressionV2KmerInfo {
    ///
    /// copy num of the node (k-mer)
    ///
    copy_num: CopyNum,
    ///
    /// frequency of the node
    ///
    freq: Freq,
    ///
    /// frequency of the intersection
    ///
    freq_intersection: Freq,
    ///
    /// total frequency of Begin->node
    ///
    freq_init: Freq,
    ///
    /// current genome size
    ///
    copy_num_total: CopyNum,
    ///
    /// current size of the intersection
    ///
    copy_num_intersection: CopyNum,
    ///
    /// expected genome size
    ///
    copy_num_total_expected: CopyNum,
    ///
    /// lambda
    ///
    penalty_weight: f64,
}

impl CompressionV2KmerInfo {
    fn new(
        copy_num: CopyNum,
        freq: Freq,
        freq_intersection: Freq,
        freq_init: Freq,
        copy_num_total: CopyNum,
        copy_num_intersection: CopyNum,
        copy_num_total_expected: CopyNum,
        penalty_weight: f64,
    ) -> Self {
        CompressionV2KmerInfo {
            copy_num,
            freq,
            freq_intersection,
            freq_init,
            copy_num_total,
            copy_num_intersection,
            copy_num_total_expected,
            penalty_weight,
        }
    }
}

pub type KmerInfos = NodeVec<DenseStorage<CompressionV2KmerInfo>>;

pub type SimpleEDbgEdgeWithKmerInfos<K> = SimpleEDbgEdgeWithAttr<K, CompressionV2KmerInfo>;

impl<K: KmerLike> FlowEdge for SimpleEDbgEdgeWithKmerInfos<K> {
    fn demand(&self) -> usize {
        0
    }
    fn capacity(&self) -> usize {
        MAX_COPY_NUM_OF_EDGE
    }
}

impl<K: KmerLike> ConvexCost for SimpleEDbgEdgeWithKmerInfos<K> {
    fn convex_cost(&self, flow: usize) -> f64 {
        unimplemented!();
    }
}

///
/// Construct NodeVec of V2KmerInfo
/// from necessary informations
///
fn create_kmer_infos<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
) -> KmerInfos {
    let mut ki = KmerInfos::new(dbg.n_nodes(), CompressionV2KmerInfo::default());

    let copy_num_total = dbg.genome_size();
    let freq_init = init_freqs.sum();

    for intersection in dbg.iter_flow_intersections(edge_freqs) {
        let copy_num_intersection = intersection.total_in_copy_num();
        let freq_intersection = intersection.total_freq().unwrap();

        // any node belongs to an intersection
        // so these loops should enumerate all nodes.
        for node_info in intersection.iter_out_nodes() {
            let node = node_info.index;
            let copy_num = node_info.copy_num;

            let freq_parents: f64 = dbg.parents(node).map(|(e, _, _)| edge_freqs[e]).sum();
            let freq = init_freqs[node] + freq_parents;

            ki[node] = CompressionV2KmerInfo::new(
                copy_num,
                freq,
                freq_intersection,
                freq_init,
                copy_num_total,
                copy_num_intersection,
                genome_size,
                penalty_weight,
            );
        }
    }

    ki
}

///
/// E-step of compression_v2
///
/// calculate edge_freqs (freq between v->w) and init_freqs (freq between Begin->w)
///
fn e_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
) -> (EdgeFreqs, NodeFreqs, Prob) {
    let phmm = dbg.to_phmm(params.clone());
    phmm.to_edge_and_init_freqs_parallel(reads)
}

///
/// M-step of compression_v2
///
/// * Construct edbg whose edge has CompressionV2KmerInfo
/// * Solve min-cost-flow
///
fn m_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
) -> (NodeCopyNums, Cost) {
    let infos = create_kmer_infos(dbg, edge_freqs, init_freqs, genome_size, penalty_weight);
    let edbg = dbg.to_edbg_with_attr(Some(&infos));
    // TODO starts from current copy nums
    let flow = min_cost_flow_convex_fast(&edbg.graph);
    match flow {
        None => panic!("compression::v2::m_step cannot find optimal flow."),
        // an edge in edbg corresponds to a node in dbg
        // so edgevec for edbg can be converted to nodevec for dbg.
        Some(copy_nums) => {
            let cost = total_cost(&edbg.graph, &copy_nums);
            (copy_nums.switch_index(), cost)
        }
    }
}
