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
use crate::em::compression::v2::CompressionV2KmerInfo;
use crate::min_flow::convex::ConvexCost;
use crate::min_flow::flow::FlowEdge;

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
    pub fn score(&self, copy_num: CopyNum) -> f64 {
        let diff = if copy_num == self.copy_num() {
            Diff::Nop
        } else if copy_num > self.copy_num() {
            Diff::Inc
        } else {
            Diff::Dec
        };
        self.x(diff) + self.y(diff) + self.z(diff) + self.r(diff)
    }
    ///
    /// `X_l log (copy_num)` term
    ///
    pub fn x(&self, diff: Diff) -> f64 {
        let new_copy_num = match diff {
            Diff::Inc => self.0.copy_num + 1,
            Diff::Dec => self.0.copy_num - 1,
            Diff::Nop => self.0.copy_num,
        };
        self.0.freq * clamped_log(new_copy_num)
    }
    ///
    /// `Y_i log (copy_num_intersection)` term
    ///
    pub fn y(&self, diff: Diff) -> f64 {
        let new_copy_num_intersection = match diff {
            Diff::Inc => self.0.copy_num_intersection + 1,
            Diff::Dec => self.0.copy_num_intersection - 1,
            Diff::Nop => self.0.copy_num_intersection,
        };
        -self.0.freq_intersection * clamped_log(new_copy_num_intersection)
    }
    ///
    /// `Z log (copy_num_total)` term
    ///
    pub fn z(&self, diff: Diff) -> f64 {
        let new_copy_num_total = match diff {
            Diff::Inc => self.0.copy_num_total + 1,
            Diff::Dec => self.0.copy_num_total - 1,
            Diff::Nop => self.0.copy_num_total,
        };
        -self.0.freq_init * clamped_log(new_copy_num_total)
    }
    ///
    /// regularization term `R`
    ///
    pub fn r(&self, diff: Diff) -> f64 {
        let new_copy_num_total = match diff {
            Diff::Inc => self.0.copy_num_total + 1,
            Diff::Dec => self.0.copy_num_total - 1,
            Diff::Nop => self.0.copy_num_total,
        };
        -self.0.penalty_weight * (new_copy_num_total - self.0.copy_num_total_expected).pow(2) as f64
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
        self.attribute.score(flow)
    }
}

// ** m-step **
// Convert
//
// ** e-step **
// Convert to min-flow network and solve it to find the best improvement of copy nums.
//
