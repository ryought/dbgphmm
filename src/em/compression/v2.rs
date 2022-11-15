//!
//! CompressionV2
//!
pub use super::kmer_info::{create_kmer_infos as create_plain_kmer_infos, KmerInfo};
use super::q::q_score;
use crate::common::{CopyNum, Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::dbg::draft::MAX_COPY_NUM_OF_EDGE;
use crate::dbg::edge_centric::impls::SimpleEDbgEdgeWithAttr;
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
pub struct CompressionV2KmerInfo(KmerInfo);

pub type V2KmerInfos = NodeVec<DenseStorage<CompressionV2KmerInfo>>;
pub type SimpleEDbgEdgeWithV2KmerInfos<K> = SimpleEDbgEdgeWithAttr<K, CompressionV2KmerInfo>;

impl std::fmt::Display for CompressionV2KmerInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{} x={:.4} y={:.4} z={:.4} 0={} +1={} -1={}",
            self.0,
            self.x(),
            self.y(),
            self.z(),
            self.score(self.0.copy_num),
            if self.0.copy_num < MAX_COPY_NUM_OF_EDGE {
                format!(
                    "{}",
                    self.score(self.0.copy_num + 1) - self.score(self.0.copy_num)
                )
            } else {
                "x".to_string()
            },
            if self.0.copy_num > 0 {
                format!(
                    "{}",
                    self.score(self.0.copy_num - 1) - self.score(self.0.copy_num)
                )
            } else {
                "x".to_string()
            },
        )
    }
}

impl CompressionV2KmerInfo {
    /// Coefficient of log(copy_num)
    pub fn x(&self) -> f64 {
        self.0.freq
    }
    /// Coefficient of (copy_num)^2
    pub fn y(&self) -> f64 {
        self.0.penalty_weight * self.0.copy_num_total as f64 / self.0.copy_num as f64
    }
    /// Coefficient of (copy_num)
    pub fn z(&self) -> f64 {
        self.0.freq_intersection / self.0.copy_num_intersection as f64
            + self.0.freq_init / self.0.copy_num_total as f64
            + self.0.penalty_weight * self.0.copy_num_total_expected as f64
    }
    ///
    ///
    pub fn score(&self, copy_num: CopyNum) -> f64 {
        if self.0.is_emittable {
            let x = self.x();
            let y = self.y();
            let z = self.z();
            assert!(x >= 0.0);
            assert!(y >= 0.0);
            assert!(z >= 0.0);
            assert!(!x.is_nan());
            assert!(!y.is_nan());
            assert!(!z.is_nan());
            if copy_num == 0 {
                // to avoid Inf*0
                -x * clamped_log(copy_num)
            } else {
                (-x * clamped_log(copy_num)) + (y * copy_num.pow(2) as f64) + (z * copy_num as f64)
            }
        } else {
            0.0
        }
    }
}

impl<K: KmerLike> FlowEdge<usize> for SimpleEDbgEdgeWithV2KmerInfos<K> {
    fn demand(&self) -> usize {
        0
    }
    fn capacity(&self) -> usize {
        MAX_COPY_NUM_OF_EDGE
    }
}

impl<K: KmerLike> ConvexCost<usize> for SimpleEDbgEdgeWithV2KmerInfos<K> {
    ///
    /// convex cost for exact compression EM algorithm
    ///
    fn convex_cost(&self, flow: usize) -> f64 {
        self.attribute.score(flow)
    }
}

///
/// create plain KmerInfos and convert it to V2KmerInfos
///
fn create_kmer_infos<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
) -> V2KmerInfos {
    let mut ki = V2KmerInfos::new(dbg.n_nodes(), CompressionV2KmerInfo::default());
    let ki0 = create_plain_kmer_infos(dbg, edge_freqs, init_freqs, genome_size, penalty_weight);
    for (node, _) in dbg.nodes() {
        ki[node] = CompressionV2KmerInfo(ki0[node]);
    }
    ki
}

///
/// E-step of compression_v2
///
/// calculate edge_freqs (freq between v->w) and init_freqs (freq between Begin->w)
///
pub fn e_step<N: DbgNode, E: DbgEdge>(
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
fn m_step_once<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
) -> (NodeCopyNums, Cost) {
    // construct edbg with KmerInfo
    let infos = create_kmer_infos(dbg, edge_freqs, init_freqs, genome_size, penalty_weight);
    let edbg = dbg.to_edbg_with_attr(Some(&infos));
    // dbg.draw_with_vecs(&[&infos], &[]);
    // dbg.draw_plain_with_vecs(&[&infos], &[]);

    // min-flow optimization starts from current copy nums
    let original_copy_nums = dbg.to_node_copy_nums().switch_index();
    let copy_nums = min_cost_flow_from_convex(&edbg.graph, &original_copy_nums);
    println!("old_copy_nums={}", original_copy_nums);
    println!("new_copy_nums={}", copy_nums);
    println!("cost_old={}", total_cost(&edbg.graph, &original_copy_nums));
    println!("cost_new={}", total_cost(&edbg.graph, &copy_nums));
    let cost_diff =
        total_cost(&edbg.graph, &copy_nums) - total_cost(&edbg.graph, &original_copy_nums);
    println!("cost_diff={}", cost_diff);
    (copy_nums.switch_index(), cost_diff)
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
) -> Dbg<N, E> {
    let mut dbg = dbg.clone();
    loop {
        let (new_copy_nums, cost_diff) =
            m_step_once(&dbg, edge_freqs, init_freqs, genome_size, penalty_weight);
        if cost_diff.is_nan() || cost_diff >= 0.0 {
            // not improved
            break;
        }
        dbg.set_node_copy_nums(&new_copy_nums);
    }
    dbg
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
        let params = PHMMParams::default();
        println!("{}", dbg);
        println!("{}", dbg.n_traverse_choices());
        let (ef, nf, p) = e_step(&dbg, &reads, &params);
        println!("{}", ef);
        println!("{}", nf);

        let lambda = 1.0;
        let genome_size = 100;

        let infos = create_kmer_infos(&dbg, &ef, &nf, genome_size, lambda);
        // dbg.draw_with_vecs(&[&nf], &[&ef]);
        dbg.draw_with_vecs(&[&infos], &[]);
        // println!("{}", infos);
        let qs = q_score(&dbg, &ef, &nf, genome_size, lambda);
        println!("{:?}", qs);

        let dbg = m_step(&dbg, &ef, &nf, genome_size, lambda);
        println!("{}", dbg);

        /*
        let mut new_dbg = dbg.clone();
        let is_updated = new_dbg.set_node_copy_nums(&ncn);
        let qs = q_score(&new_dbg, &ef, &nf, genome_size, lambda);
        println!("{:?}", qs);
        */
    }
}
