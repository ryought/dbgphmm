//!
//!
//!
use crate::common::{CopyNum, Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::dbg::edge_centric::impls::{SimpleEDbgEdgeWithAttr, MAX_COPY_NUM_OF_EDGE};
use crate::hmmv2::{EdgeFreqs, NodeFreqs};
use crate::min_flow::utils::clamped_log;
use crate::vector::{DenseStorage, EdgeVec, NodeVec, Storage};

#[derive(Clone, Debug, Copy, PartialEq, Default)]
pub struct KmerInfo {
    ///
    /// This is emittable kmer or not.
    /// If not emittable, the kmer will be excluded for the cost.
    ///
    pub is_emittable: bool,
    ///
    /// copy num of the node (k-mer)
    ///
    pub copy_num: CopyNum,
    ///
    /// frequency of the node
    ///
    pub freq: Freq,
    ///
    /// frequency of the intersection
    ///
    pub freq_intersection: Freq,
    ///
    /// total frequency from Begin
    ///
    pub freq_init: Freq,
    ///
    /// current genome size
    ///
    pub copy_num_total: CopyNum,
    ///
    /// current size of the intersection
    ///
    pub copy_num_intersection: CopyNum,
    ///
    /// expected genome size
    ///
    pub copy_num_total_expected: CopyNum,
    ///
    /// lambda
    ///
    pub penalty_weight: f64,
}

impl std::fmt::Display for KmerInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "e={} c={} ci={} cG={} cG0={} f={} fi={} fB={} w={}",
            self.is_emittable,
            self.copy_num,
            self.copy_num_intersection,
            self.copy_num_total,
            self.copy_num_total_expected,
            self.freq,
            self.freq_intersection,
            self.freq_init,
            self.penalty_weight,
        )
    }
}

impl KmerInfo {
    pub fn new(
        is_emittable: bool,
        copy_num: CopyNum,
        freq: Freq,
        freq_intersection: Freq,
        freq_init: Freq,
        copy_num_total: CopyNum,
        copy_num_intersection: CopyNum,
        copy_num_total_expected: CopyNum,
        penalty_weight: f64,
    ) -> Self {
        assert!(copy_num >= 0);
        assert!(freq >= 0.0);
        assert!(freq_intersection >= 0.0);
        assert!(freq_init >= 0.0);
        assert!(copy_num_total >= 0);
        assert!(copy_num_intersection >= 0);
        assert!(copy_num_total_expected >= 0);
        assert!(penalty_weight >= 0.0);
        KmerInfo {
            is_emittable,
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

///
/// NodeVec of KmerInfo
/// can be constructed by `create_kmer_infos`.
///
pub type KmerInfos = NodeVec<DenseStorage<KmerInfo>>;

///
/// Construct NodeVec of KmerInfo from necessary informations
///
pub fn create_kmer_infos<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
) -> KmerInfos {
    let mut ki = KmerInfos::new(dbg.n_nodes(), KmerInfo::default());

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
            let is_emittable = dbg.node(node).is_emittable();

            let freq_parents: f64 = dbg.parents(node).map(|(e, _, _)| edge_freqs[e]).sum();
            let freq = init_freqs[node] + freq_parents;

            ki[node] = KmerInfo::new(
                is_emittable,
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
