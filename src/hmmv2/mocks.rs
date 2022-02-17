//!
//! Mock PHMMs for testing
//!
use super::common::PModel;
use crate::graph::mocks::*;
use crate::graph::seq_graph::SeqGraph;
use crate::hmmv2::params::PHMMParams;

///
/// Sample linear genome (10bp) "ATTCGATCGT"
///
pub fn mock_linear_phmm(param: PHMMParams) -> PModel {
    mock_linear().to_seq_graph().to_phmm(param)
}

///
/// Create random linear phmm
///
pub fn mock_linear_random_phmm(length: usize, seed: u64, param: PHMMParams) -> PModel {
    mock_linear_random(length, seed)
        .to_seq_graph()
        .to_phmm(param)
}

/// Example small grpah
///
/// ```text
///     a(x2) ----   ------ c(x2)
///               \ /
///                X
///               / \
///     b(x2) ----   ------ d(x2)
/// ```
///
pub fn mock_crossing_phmm(has_edge_copy_number: bool, param: PHMMParams) -> PModel {
    mock_crossing(has_edge_copy_number)
        .to_seq_graph()
        .to_phmm(param)
}
