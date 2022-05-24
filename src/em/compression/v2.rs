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
    freq_node: Freq,
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
    size_total: CopyNum,
    ///
    /// current size of the intersection
    ///
    size_intersection: CopyNum,
    ///
    /// expected genome size
    ///
    size_total_expected: CopyNum,
    ///
    /// lambda
    ///
    penalty_weight: f64,
}
