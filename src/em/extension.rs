//!
//! Extension
//!
//! ## E-step
//!
//!
//! ## M-step
//!
//!
use crate::common::{CopyNum, Freq};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, EdgeCopyNums};
use crate::graph::Bipartite;
use crate::hmmv2::freq::Reads;
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::trans_table::EdgeFreqs;
use crate::kmer::kmer::KmerLike;
use crate::min_flow::convex::ConvexCost;
use crate::min_flow::flow::FlowEdge;
use crate::min_flow::utils::clamped_log;
use petgraph::dot::Dot;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

///
/// Extension algorithm
///
/// ## Details
///
/// * e-step
///     Calculate edge_freqs on dbg.
///
/// * m-step
///     Maximize the score for each intersections.
///
/// ## TODOs
///
/// * avoid dbg copy?
///
pub fn extension<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
) -> Dbg<N, E> {
    unimplemented!();
}

///
/// E-step of extension
///
/// ## Details
///
/// * convert dbg into phmm.
/// * calculate edge frequencies by forward/backward algorithm to emit the reads.
///
fn extension_e_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
) -> EdgeFreqs {
    let phmm = dbg.to_phmm(params.clone());
    phmm.to_edge_freqs(reads)
}

///
/// M-step of extension
///
/// ## Params
///
/// * `dbg`
///     de bruijn graph
/// * `edge_freqs`
///     edge usage frequencies.
///
/// ## Details
///
///
fn extension_m_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
) -> EdgeCopyNums {
    unimplemented!();
}

//
// Definition of Augmented-intersection information collection
//

/// Node info
#[derive(Clone, Debug)]
pub struct FlowIntersectionNode {
    index: NodeIndex,
    copy_num: CopyNum,
}

impl FlowIntersectionNode {
    /// Constructor
    pub fn new(index: NodeIndex, copy_num: CopyNum) -> Self {
        FlowIntersectionNode { index, copy_num }
    }
}

impl std::fmt::Display for FlowIntersectionNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}(x{})", self.index.index(), self.copy_num)
    }
}

/// Edge between in-node and out-node
#[derive(Clone, Debug)]
pub struct FlowIntersectionEdge {
    index: EdgeIndex,
    freq: Freq,
}

impl FlowIntersectionEdge {
    /// Constructor
    pub fn new(index: EdgeIndex, freq: Freq) -> Self {
        FlowIntersectionEdge { index, freq }
    }
}

impl std::fmt::Display for FlowIntersectionEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}({})", self.index.index(), self.freq)
    }
}

#[derive(Clone, Debug)]
pub struct FlowIntersection<K: KmerLike> {
    bi: Bipartite<K, FlowIntersectionNode, FlowIntersectionEdge>,
}

impl<K: KmerLike> FlowIntersection<K> {
    pub fn new(
        km1mer: K,
        in_nodes: Vec<FlowIntersectionNode>,
        out_nodes: Vec<FlowIntersectionNode>,
        edges: Vec<FlowIntersectionEdge>,
    ) -> Self {
        FlowIntersection {
            bi: Bipartite::new(km1mer, in_nodes, out_nodes, edges),
        }
    }
    pub fn from<F: Fn(usize, usize) -> FlowIntersectionEdge>(
        km1mer: K,
        in_nodes: Vec<FlowIntersectionNode>,
        out_nodes: Vec<FlowIntersectionNode>,
        edge_fn: F,
    ) -> Self {
        FlowIntersection {
            bi: Bipartite::from(km1mer, in_nodes, out_nodes, edge_fn),
        }
    }
}

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
    Spanning { freq: Freq, edge_index: EdgeIndex },
}

impl FlowEdge for IntersectionGraphEdge {
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

impl ConvexCost for IntersectionGraphEdge {
    fn convex_cost(&self, flow: usize) -> f64 {
        match self {
            IntersectionGraphEdge::Spanning { freq, .. } => -(*freq) * clamped_log(flow),
            _ => 0.0,
        }
    }
}

impl<K: KmerLike> FlowIntersection<K> {
    /// Get optimized copy numbers of edges.
    /// by converting the bipartite into flow network definitions
    ///
    /// ## TODOs
    /// * what is the return type?
    pub fn optimize(&self) {
        unimplemented!();
    }
    ///
    /// Convert FlowIntersection into DiGraph
    /// that can be solved as a convex min flow problem.
    ///
    pub fn to_flow_graph(&self) -> IntersectionGraph {
        let mut g = IntersectionGraph::new();

        // (1) Node
        // (1-a) source and terminal
        let s = g.add_node(IntersectionGraphNode::Source);
        let t = g.add_node(IntersectionGraphNode::Terminal);
        // (1-b) left/right nodes
        let vs: Vec<NodeIndex> = self
            .bi
            .iter_in_nodes()
            .map(|v| g.add_node(IntersectionGraphNode::Left(v.index)))
            .collect();
        let ws: Vec<NodeIndex> = self
            .bi
            .iter_out_nodes()
            .map(|w| g.add_node(IntersectionGraphNode::Right(w.index)))
            .collect();

        // (2) Edge
        // (2-a) LoopBack
        g.add_edge(t, s, IntersectionGraphEdge::LoopBack);
        // (2-b) Node
        for (i, n) in self.bi.iter_in_nodes().enumerate() {
            g.add_edge(
                s,
                vs[i],
                IntersectionGraphEdge::Node {
                    copy_num: n.copy_num,
                },
            );
        }
        for (i, n) in self.bi.iter_out_nodes().enumerate() {
            g.add_edge(
                ws[i],
                t,
                IntersectionGraphEdge::Node {
                    copy_num: n.copy_num,
                },
            );
        }
        // (2-c) Spanning
        for i in 0..self.bi.n_in() {
            for j in 0..self.bi.n_out() {
                let e = self.bi.edge(i, j);
                g.add_edge(
                    vs[i],
                    ws[j],
                    IntersectionGraphEdge::Spanning {
                        freq: e.freq,
                        edge_index: e.index,
                    },
                );
            }
        }

        g
    }
}

impl<K: KmerLike> std::fmt::Display for FlowIntersection<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.bi)
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};
    use crate::kmer::veckmer::VecKmer;

    #[test]
    fn flow_intersection_construction() {
        let in_nodes = vec![
            FlowIntersectionNode::new(ni(0), 2),
            FlowIntersectionNode::new(ni(1), 4),
        ];
        let out_nodes = vec![
            FlowIntersectionNode::new(ni(3), 1),
            FlowIntersectionNode::new(ni(4), 5),
        ];
        let edges = vec![
            FlowIntersectionEdge::new(ei(0), 5.1),
            FlowIntersectionEdge::new(ei(1), 5.1),
            FlowIntersectionEdge::new(ei(2), 5.1),
            FlowIntersectionEdge::new(ei(3), 5.1),
        ];
        let fi = FlowIntersection::new(VecKmer::from_bases(b"TCG"), in_nodes, out_nodes, edges);
        println!("{}", fi);
        let g = fi.to_flow_graph();
        println!("{:?}", Dot::with_config(&g, &[]));
    }
}
