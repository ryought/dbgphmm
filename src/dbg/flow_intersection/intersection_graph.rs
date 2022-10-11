use crate::common::{CopyNum, Freq};
use crate::min_flow::convex::ConvexCost;
use crate::min_flow::flow::FlowEdge;
use crate::min_flow::utils::clamped_log;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

const MAX_COPY_NUM_OF_EDGE: usize = 1000;
pub type IntersectionGraph = DiGraph<IntersectionGraphNode, IntersectionGraphEdge>;

#[derive(Clone, Debug)]
pub enum IntersectionGraphNode {
    Source,
    Terminal,
    Left(NodeIndex),
    Right(NodeIndex),
}

#[derive(Clone, Debug)]
pub enum IntersectionGraphEdge {
    /// Loop back edge from t (terminal) to s (source).
    ///
    /// demand=0, capacity=inf, cost=0
    LoopBack,
    /// Node edge
    /// * from s to v (left nodes in bipartite)
    /// * from w (right nodes in bipartite) to t
    ///
    /// demand=copy_num, capacity=copy_num, cost=0
    Node { copy_num: CopyNum },
    /// Spanning edge between v and w (left/right nodes in bipartite)
    ///
    /// demand=0, capacity=inf, cost=-freq*log(flow)
    Spanning {
        freq: Freq,
        edge_index: EdgeIndex,
        edge_index_in_intersection: usize,
    },
}

impl FlowEdge<usize> for IntersectionGraphEdge {
    fn demand(&self) -> usize {
        match self {
            IntersectionGraphEdge::Node { copy_num } => *copy_num,
            _ => 0,
        }
    }
    fn capacity(&self) -> usize {
        match self {
            IntersectionGraphEdge::Node { copy_num } => *copy_num,
            _ => MAX_COPY_NUM_OF_EDGE,
        }
    }
}

impl ConvexCost<usize> for IntersectionGraphEdge {
    fn convex_cost(&self, flow: usize) -> f64 {
        match self {
            IntersectionGraphEdge::Spanning { freq, .. } => -(*freq) * clamped_log(flow),
            _ => 0.0,
        }
    }
}
