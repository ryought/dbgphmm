//!
//!
//!
use crate::common::{CopyNum, Freq, Reads};
use crate::dbg::edge_centric::impls::{SimpleEDbgEdgeWithAttr, MAX_COPY_NUM_OF_EDGE};
use crate::min_flow::utils::clamped_log;

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
    /// total frequency from Begin
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
            "e={} c={} ci={} cG={} cG0={} f={} fi={} fB={} w={} x={:.4} y={:.4} z={:.4} 0={} +1={} -1={}",
            self.is_emittable,
            self.copy_num,
            self.copy_num_intersection,
            self.copy_num_total,
            self.copy_num_total_expected,
            self.freq,
            self.freq_intersection,
            self.freq_init,
            self.penalty_weight,
            self.x(),
            self.y(),
            self.z(),
            self.score(self.copy_num),
            if self.copy_num < MAX_COPY_NUM_OF_EDGE {
                format!(
                    "{}",
                    self.score(self.copy_num + 1) - self.score(self.copy_num)
                )
            } else {
                "x".to_string()
            },
            if self.copy_num > 0 {
                format!(
                    "{}",
                    self.score(self.copy_num - 1) - self.score(self.copy_num)
                )
            } else {
                "x".to_string()
            },
        )
    }
}

impl CompressionV2KmerInfo {
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
