//!
//! CompressionV2
//!
use crate::common::{CopyNum, Freq};

#[derive(Clone, Debug)]
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
