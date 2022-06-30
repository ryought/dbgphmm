//!
//! Compression V3
//! By Naive main factor improving
//!
use super::kmer_info::{create_kmer_infos as create_plain_kmer_infos, KmerInfo};
use super::q::{q_score_clamped, QScore};
use crate::common::{CopyNum, Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::{EdgeFreqs, NodeFreqs};
use crate::kmer::kmer::KmerLike;
use crate::min_flow::residue::improve_flow_convex;
use crate::min_flow::utils::clamped_log;
use crate::min_flow::{min_cost_flow_from_convex, total_cost, Cost};
use crate::prob::Prob;
use crate::vector::{DenseStorage, EdgeVec, NodeVec, Storage};

use crate::dbg::edge_centric::impls::{SimpleEDbgEdgeWithAttr, MAX_COPY_NUM_OF_EDGE};
use crate::min_flow::convex::ConvexCost;
use crate::min_flow::flow::FlowEdge;

// use same e-step function in V2
pub use crate::em::compression::v2::e_step;

#[derive(Clone, Debug, Copy, PartialEq, Default)]
pub struct CompressionV3KmerInfo(KmerInfo);

pub type V3KmerInfos = NodeVec<DenseStorage<CompressionV3KmerInfo>>;
pub type SimpleEDbgEdgeWithV3KmerInfos<K> = SimpleEDbgEdgeWithAttr<K, CompressionV3KmerInfo>;

#[derive(Clone, Debug, Copy, PartialEq)]
pub enum Diff {
    /// +1
    Inc,
    /// -1
    Dec,
    /// 0
    Nop,
}

impl Diff {
    ///
    /// Apply the diff operation
    ///
    pub fn apply(&self, value: usize) -> usize {
        match self {
            Diff::Inc => value.saturating_add(1),
            Diff::Dec => value.saturating_sub(1),
            Diff::Nop => value,
        }
    }
}

impl CompressionV3KmerInfo {
    ///
    /// Current copy number of the kmer
    ///
    pub fn copy_num(&self) -> CopyNum {
        self.0.copy_num
    }
    ///
    /// Exact Q function score
    ///
    pub fn score_from_copy_num(&self, copy_num: CopyNum) -> f64 {
        let diff = if copy_num == self.copy_num() {
            Diff::Nop
        } else if copy_num > self.copy_num() {
            Diff::Inc
        } else {
            Diff::Dec
        };
        self.score(diff)
    }
    ///
    /// Exact Q function score
    ///
    pub fn score(&self, diff: Diff) -> f64 {
        self.x(diff) + self.y(diff) + self.z(diff) + self.r(diff)
    }
    ///
    /// `X_l log (copy_num)` term
    ///
    pub fn x(&self, diff: Diff) -> f64 {
        let new_copy_num = diff.apply(self.0.copy_num);
        -self.0.freq * clamped_log(new_copy_num)
    }
    ///
    /// `Y_i log (copy_num_intersection)` term
    ///
    pub fn y(&self, diff: Diff) -> f64 {
        let new_copy_num_intersection = diff.apply(self.0.copy_num_intersection);
        self.0.freq_intersection * clamped_log(new_copy_num_intersection)
    }
    ///
    /// `Z log (copy_num_total)` term
    ///
    pub fn z(&self, diff: Diff) -> f64 {
        let new_copy_num_total = diff.apply(self.0.copy_num_total);
        self.0.freq_init * clamped_log(new_copy_num_total)
    }
    ///
    /// regularization term `R`
    ///
    pub fn r(&self, diff: Diff) -> f64 {
        let new_copy_num_total = diff.apply(self.0.copy_num_total);
        self.0.penalty_weight
            * new_copy_num_total
                .abs_diff(self.0.copy_num_total_expected)
                .pow(2) as f64
    }
}

impl std::convert::From<KmerInfo> for CompressionV3KmerInfo {
    fn from(kmer_info: KmerInfo) -> Self {
        CompressionV3KmerInfo(kmer_info)
    }
}

impl std::fmt::Display for CompressionV3KmerInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "+1={:.4}({:.4},{:.4},{:.4},{:.4}) -1={:.4}({:.4},{:.4},{:.4},{:.4}) ({})",
            // +1
            self.score(Diff::Inc) - self.score(Diff::Nop),
            self.x(Diff::Inc) - self.x(Diff::Nop),
            self.y(Diff::Inc) - self.y(Diff::Nop),
            self.z(Diff::Inc) - self.z(Diff::Nop),
            self.r(Diff::Inc) - self.r(Diff::Nop),
            // -1
            self.score(Diff::Dec) - self.score(Diff::Nop),
            self.x(Diff::Dec) - self.x(Diff::Nop),
            self.y(Diff::Dec) - self.y(Diff::Nop),
            self.z(Diff::Dec) - self.z(Diff::Nop),
            self.r(Diff::Dec) - self.r(Diff::Nop),
            self.0,
        )
    }
}

impl<K: KmerLike> FlowEdge for SimpleEDbgEdgeWithV3KmerInfos<K> {
    ///
    /// demand is set to be current copy_num - 1
    ///
    fn demand(&self) -> usize {
        0_usize.max(self.attribute.copy_num() - 1)
    }
    ///
    /// capacity is set to be current copy_num + 1
    ///
    fn capacity(&self) -> usize {
        MAX_COPY_NUM_OF_EDGE.min(self.attribute.copy_num() + 1)
    }
}

impl<K: KmerLike> ConvexCost for SimpleEDbgEdgeWithV3KmerInfos<K> {
    ///
    ///
    ///
    fn convex_cost(&self, flow: usize) -> f64 {
        self.attribute.score_from_copy_num(flow)
    }
}

///
/// create plain KmerInfos and convert it to V3KmerInfos
///
fn create_kmer_infos<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
) -> V3KmerInfos {
    let mut ki = V3KmerInfos::new(dbg.n_nodes(), CompressionV3KmerInfo::default());
    let ki0 = create_plain_kmer_infos(dbg, edge_freqs, init_freqs, genome_size, penalty_weight);
    for (node, _) in dbg.nodes() {
        ki[node] = CompressionV3KmerInfo::from(ki0[node]);
    }
    ki
}

///
/// M-step of compression_v3
///
/// Returned value is (new copy_nums, expected improvement of the cost)
///
/// find a candidate update of current flow
/// using variational approximation
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
    // dbg.draw_with_vecs(&[&infos], &[]);
    // dbg.draw_plain_with_vecs(&[&infos], &[]);

    // min-flow optimization starts from current copy nums
    let original_copy_nums = dbg.to_node_copy_nums().switch_index();
    let (copy_nums, cost_diff) = match improve_flow_convex(&edbg.graph, &original_copy_nums) {
        Some(copy_nums) => {
            let cost_diff =
                total_cost(&edbg.graph, &copy_nums) - total_cost(&edbg.graph, &original_copy_nums);
            (copy_nums, cost_diff)
        }
        None => (original_copy_nums, 0.0),
    };

    (copy_nums.switch_index(), cost_diff)
}

pub fn compression_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    genome_size: CopyNum,
    penalty_weight: f64,
) -> (Dbg<N, E>, bool, CompressionV3Log<N, E>) {
    // [1] e-step
    let (edge_freqs, init_freqs, p) = e_step(dbg, reads, params);

    // calculate current q-score
    let q_score = q_score_clamped(dbg, &edge_freqs, &init_freqs, genome_size, penalty_weight);

    // [2] m-step
    let (copy_nums_new, cost_diff) =
        m_step(dbg, &edge_freqs, &init_freqs, genome_size, penalty_weight);

    // construct new candidate dbg
    let mut new_dbg = dbg.clone();
    let is_updated = new_dbg.set_node_copy_nums(&copy_nums_new);
    // new_dbg.remove_zero_copy_node();
    let q_score_new = q_score_clamped(
        &new_dbg,
        &edge_freqs,
        &init_freqs,
        genome_size,
        penalty_weight,
    );

    // [3] history
    let log = CompressionV3Log::new(q_score, q_score_new, cost_diff, is_updated, new_dbg.clone());
    println!("{}", log);

    (new_dbg, is_updated, log)
}

///
/// Log information store of each iteration in compression v3
///
#[derive(Clone)]
pub struct CompressionV3Log<N: DbgNode, E: DbgEdge> {
    /// q-score before m-step
    pub q0: QScore,
    /// q-score after m-step
    pub q1: QScore,
    /// cost improvement with variational-approximated q-score
    pub cost_diff: Cost,
    /// the dbg changed by this compression step?
    pub is_updated: bool,
    /// resulting dbg
    pub dbg: Dbg<N, E>,
}

impl<N: DbgNode, E: DbgEdge> CompressionV3Log<N, E> {
    pub fn new(q0: QScore, q1: QScore, cost_diff: Cost, is_updated: bool, dbg: Dbg<N, E>) -> Self {
        CompressionV3Log {
            q0,
            q1,
            cost_diff,
            is_updated,
            dbg,
        }
    }
}

impl<N: DbgNode, E: DbgEdge> std::fmt::Display for CompressionV3Log<N, E> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\tis_updated={}\t{}",
            self.q0, self.q1, self.cost_diff, self.is_updated, self.dbg
        )
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks::*;

    #[test]
    fn em_compression_v3_m() {
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

        let lambda = 0.0;
        let genome_size = 9;

        let infos = create_kmer_infos(&dbg, &ef, &nf, genome_size, lambda);
        // dbg.draw_with_vecs(&[&nf], &[&ef]);
        dbg.draw_with_vecs(&[&infos], &[]);
        // println!("{}", infos);
        // let qs = q_score(&dbg, &ef, &nf, genome_size, lambda);
        // println!("{:?}", qs);

        let (copy_nums, d) = m_step(&dbg, &ef, &nf, genome_size, lambda);
        println!("{}", copy_nums);
        println!("{}", d);
        assert_eq!(
            copy_nums.to_vec(),
            vec![0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0]
        );
    }

    #[test]
    fn em_compression_v3_step() {
        let dbg = mock_intersection();
        let reads = Reads {
            reads: vec![
                b"AACTAGCTT".to_vec(),
                b"AACTAGCTT".to_vec(),
                b"AACTAGCTT".to_vec(),
            ],
        };
        let params = PHMMParams::default();
        println!("dbg0={}", dbg);

        let (new_dbg, is_updated, log) = compression_step(&dbg, &reads, &params, 9, 0.0);
        println!("dbg1={}", new_dbg);
    }
}
