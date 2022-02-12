//!
//! Mock PHMMs for testing
//!
use super::common::PModel;
use crate::graph::mocks::*;
use crate::graph::seq_graph::SeqGraph;
use crate::hmm::params::PHMMParams;

///
/// Sample linear genome (10bp) "ATTCGATCGT"
///
pub fn mock_linear_phmm(param: PHMMParams) -> PModel {
    mock_linear().to_seq_graph().to_phmm(param)
}
