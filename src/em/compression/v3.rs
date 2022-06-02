//!
//! Compression V3
//! By Naive main factor improving
//!
use super::kmer_info::{create_kmer_infos as create_plain_kmer_infos, KmerInfo};
use crate::common::{CopyNum, Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::{EdgeFreqs, NodeFreqs};
use crate::kmer::kmer::KmerLike;
use crate::min_flow::utils::clamped_log;
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
            "+1={}({},{},{},{}) -1={}({},{},{},{}) ({})",
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

// ** m-step **
// Convert
//
// ** e-step **
// Convert to min-flow network and solve it to find the best improvement of copy nums.
//

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks::*;

    #[test]
    fn em_compression_v3() {
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

        // let dbg = m_step(&dbg, &ef, &nf, genome_size, lambda);
        // println!("{}", dbg);
    }
}
