//!
//! Constructor of draft dbg
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::common::{CopyNum, Freq, Seq};
use crate::hmmv2::freq::NodeFreqs;
use crate::kmer::kmer::KmerLike;
use crate::min_flow::{
    convex::ConvexCost, flow::FlowEdge, min_cost_flow_convex_fast, total_cost, Cost,
};

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    pub fn create_draft_from_seqs<T>(k: usize, seqs: T) -> Self
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        unimplemented!();
    }
    ///
    /// by solving min-cost-flow
    ///
    pub fn min_squared_error_copy_nums_from_freqs(
        &self,
        freqs: &NodeFreqs,
    ) -> Option<(NodeCopyNums, Cost)> {
        let graph = self.to_edbg_graph(
            |_| (),
            |v, weight| {
                MinSquaredErrorCopyNumAndFreq::new(
                    weight.is_emittable(),
                    weight.copy_num(),
                    freqs[v],
                )
            },
        );
        min_cost_flow_convex_fast(&graph).map(|flow| {
            let cost = total_cost(&graph, &flow);
            // an edge in edbg corresponds to a node in dbg
            // so edgevec for edbg can be converted to nodevec for dbg.
            (flow.switch_index(), cost)
        })
    }
}

///
/// Edge attribute for min_squared_error_copy_nums_from_freqs
///
/// FlowEdge
/// * demand = 0
/// * capacity = +inf
///
/// ConvexCost
/// * cost = |c - f|^2
///
#[derive(Clone, Debug)]
struct MinSquaredErrorCopyNumAndFreq {
    is_target: bool,
    copy_num: CopyNum,
    freq: Freq,
}

impl MinSquaredErrorCopyNumAndFreq {
    ///
    /// constructor
    ///
    pub fn new(is_target: bool, copy_num: CopyNum, freq: Freq) -> Self {
        MinSquaredErrorCopyNumAndFreq {
            is_target,
            copy_num,
            freq,
        }
    }
}

///
/// maximum copy number
///
/// this corresponds to the capacity of edbg min-flow calculation.
///
pub const MAX_COPY_NUM_OF_EDGE: usize = 1000;

impl FlowEdge<usize> for MinSquaredErrorCopyNumAndFreq {
    fn demand(&self) -> usize {
        0
    }
    fn capacity(&self) -> usize {
        MAX_COPY_NUM_OF_EDGE
    }
}

///
/// Use edbg edge (with a freq) in min-flow.
///
/// if the kmer corresponding to the edge is not emittable, the cost
/// should be ignored.
///
impl ConvexCost<usize> for MinSquaredErrorCopyNumAndFreq {
    fn convex_cost(&self, copy_num: usize) -> f64 {
        if self.is_target {
            (copy_num as f64 - self.freq).powi(2)
        } else {
            0.0
        }
    }
}
