//!
//! Calculate Node/Edge (i.e. hidden state/transition) usage frequencies
//! from the result of Forward/Backward.
//!
use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::table::{PHMMResult, PHMMTable};
use crate::common::Freq;
use crate::prob::Prob;
use crate::vector::{DenseStorage, NodeVec, Storage};

/// Frequency (f64) assigned to each nodes
/// It cannot assume the sparcity, so we use the dense storage
/// TODO
pub type NodeFreq = NodeVec<DenseStorage<Freq>>;

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    /// Calculate the node usage frequency
    ///
    /// `freq[i][t_k]`
    /// = P(base `x[i]` is emitted from node `t_k`)
    /// = `ft_i[k] * bt_i+1[k]`
    ///
    /// * `t` is a type of state, either Match, Ins, Del
    /// * `k` is a node index
    ///
    pub fn to_node_freq<S>(&self, forward: &PHMMResult<S>, backward: &PHMMResult<S>) -> NodeFreq
    where
        S: Storage<Item = Prob>,
    {
        unimplemented!();
    }
    fn node_freq<S>(&self, f0: &PHMMTable<S>, b1: &PHMMTable<S>) -> PHMMTable<S>
    where
        S: Storage<Item = Prob>,
    {
        unimplemented!();
    }
    // pub fn to_edge_freq<S>() {}
}

//
// Tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::mocks::mock_linear;
    use crate::hmm::params::PHMMParams;
    use crate::vector::DenseStorage;
    #[test]
    fn hmm_freq_mock_linear() {
        let phmm = mock_linear()
            .to_seq_graph()
            .to_phmm(PHMMParams::high_error());
        let rf: PHMMResult<DenseStorage<Prob>> = phmm.forward(b"CGATC");
        let rb: PHMMResult<DenseStorage<Prob>> = phmm.backward(b"CGATC");
        phmm.to_node_freq(&rf, &rb);
    }
}
