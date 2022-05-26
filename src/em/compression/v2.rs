//!
//! CompressionV2
//!
use crate::common::{CopyNum, Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::dbg::edge_centric::impls::{SimpleEDbgEdgeWithAttr, MAX_COPY_NUM_OF_EDGE};
use crate::hmmv2::common::{PHMMEdge, PHMMNode};
use crate::hmmv2::freq::NodeFreqs;
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::trans_table::EdgeFreqs;
use crate::kmer::kmer::{Kmer, KmerLike};
use crate::min_flow::convex::ConvexCost;
use crate::min_flow::flow::FlowEdge;
use crate::min_flow::utils::clamped_log;
use crate::min_flow::{min_cost_flow_from_convex, total_cost, Cost};
use crate::prob::Prob;
use crate::vector::{DenseStorage, EdgeVec, NodeVec, Storage};

#[derive(Clone, Debug, Copy, PartialEq, Default)]
pub struct CompressionV2KmerInfo {
    ///
    /// This is emittable kmer or not.
    /// If not emittable, the kmer will be excluded for the cost.
    ///
    is_emittable: bool,
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

impl std::fmt::Display for CompressionV2KmerInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "c={} f={} fi={} fB={} cG={} ci={} cG0={} w={} x={} y={} z={} 0={} +1={} -1={}",
            self.copy_num,
            self.freq,
            self.freq_intersection,
            self.freq_init,
            self.copy_num_total,
            self.copy_num_intersection,
            self.copy_num_total_expected,
            self.penalty_weight,
            self.x(),
            self.y(),
            self.z(),
            self.score(self.copy_num),
            self.score(self.copy_num + 1) - self.score(self.copy_num),
            self.score(self.copy_num - 1) - self.score(self.copy_num),
        )
    }
}

impl CompressionV2KmerInfo {
    fn new(
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
        CompressionV2KmerInfo {
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
    /// Coefficient of log(copy_num)
    pub fn x(&self) -> f64 {
        self.freq
    }
    /// Coefficient of (copy_num)^2
    pub fn y(&self) -> f64 {
        self.penalty_weight * self.copy_num_total as f64 / self.copy_num as f64
    }
    /// Coefficient of (copy_num)
    pub fn z(&self) -> f64 {
        self.freq_intersection / self.copy_num_intersection as f64
            + self.freq_init / self.copy_num_total as f64
            + self.penalty_weight * self.copy_num_total_expected as f64
    }
    pub fn score(&self, copy_num: CopyNum) -> f64 {
        if self.is_emittable {
            let x = self.x();
            let y = self.y();
            let z = self.z();
            assert!(x >= 0.0);
            assert!(y >= 0.0);
            assert!(z >= 0.0);
            assert!(!x.is_nan());
            assert!(!y.is_nan());
            assert!(!z.is_nan());
            (-x * clamped_log(copy_num)) + (y * copy_num.pow(2) as f64) + (z * copy_num as f64)
        } else {
            0.0
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
    ///
    /// convex cost for exact compression EM algorithm
    ///
    fn convex_cost(&self, flow: usize) -> f64 {
        self.attribute.score(flow)
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
            let is_emittable = dbg.node(node).is_emittable();

            let freq_parents: f64 = dbg.parents(node).map(|(e, _, _)| edge_freqs[e]).sum();
            let freq = init_freqs[node] + freq_parents;

            ki[node] = CompressionV2KmerInfo::new(
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
    // construct edbg with KmerInfo
    let infos = create_kmer_infos(dbg, edge_freqs, init_freqs, genome_size, penalty_weight);
    let edbg = dbg.to_edbg_with_attr(Some(&infos));

    // min-flow optimization starts from current copy nums
    let original_copy_nums = dbg.to_node_copy_nums().switch_index();
    let copy_nums = min_cost_flow_from_convex(&edbg.graph, &original_copy_nums);
    println!("cost_old={}", total_cost(&edbg.graph, &original_copy_nums));
    println!("cost_new={}", total_cost(&edbg.graph, &copy_nums));
    let cost_diff =
        total_cost(&edbg.graph, &copy_nums) - total_cost(&edbg.graph, &original_copy_nums);
    (copy_nums.switch_index(), cost_diff)
}

#[derive(Clone, Debug, Copy, Default)]
struct QScore {
    init: f64,
    trans: f64,
    prior: f64,
}

impl QScore {
    pub fn total(&self) -> f64 {
        self.init + self.trans + self.prior
    }
}

///
/// Calculate (exact) Q function score.
///
fn q_score<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
) -> QScore {
    let mut qs = QScore::default();
    let phmm = dbg.to_phmm(PHMMParams::default());

    for (node, node_weight) in phmm.nodes() {
        if phmm.is_emittable(node) {
            // (1) init score
            // A(Begin, v) log p_init(v)
            qs.init += init_freqs[node] * node_weight.init_prob().to_log_value();

            for (edge, child, edge_weight) in phmm.childs(node) {
                if phmm.is_emittable(child) {
                    // (2) trans score
                    // A(v,w) log p_trans(v, w)
                    qs.trans += edge_freqs[edge] * edge_weight.trans_prob().to_log_value();
                }
            }
        }
    }

    // (3) prior score
    // -lambda (genome_size - genome_size_expected)^2
    let size_diff = genome_size as f64 - dbg.genome_size() as f64;
    qs.prior = -penalty_weight * size_diff.powi(2);

    qs
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks::*;

    #[test]
    fn em_compression_v2_e_step() {
        let dbg = mock_intersection();
        let reads = Reads {
            reads: vec![
                b"AACTAGCTT".to_vec(),
                b"AACTAGCTT".to_vec(),
                b"AACTAGCTT".to_vec(),
            ],
        };
        let params = PHMMParams::zero_error();
        println!("{}", dbg);
        println!("{}", dbg.n_traverse_choices());
        let (ef, nf, p) = e_step(&dbg, &reads, &params);
        println!("{}", ef);
        println!("{}", nf);

        let lambda = 0.0;
        let genome_size = 10;

        let infos = create_kmer_infos(&dbg, &ef, &nf, genome_size, lambda);
        dbg.draw_with_vecs(&[&nf], &[&ef]);
        // println!("{}", infos);
        let qs = q_score(&dbg, &ef, &nf, genome_size, lambda);
        println!("{:?}", qs);

        let (ncn, c) = m_step(&dbg, &ef, &nf, genome_size, lambda);
        println!("{}", ncn);

        let mut new_dbg = dbg.clone();
        let is_updated = new_dbg.set_node_copy_nums(&ncn);
        let qs = q_score(&new_dbg, &ef, &nf, genome_size, lambda);
        println!("{:?}", qs);
    }
}
