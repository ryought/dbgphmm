//!
//!
//!
use super::impls::{SimpleEDbgEdge, SimpleEDbgNode};
use super::{EDbgEdge, EDbgEdgeBase, EDbgNode};
use crate::common::{CopyNum, Freq, Sequence};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase, EdgeCopyNums, NodeCopyNums};
use crate::graph::compact::compact_simple_paths;
use crate::kmer::kmer::KmerLike;
use crate::min_flow::flow::{Flow, FlowEdge};
use crate::min_flow::{Cost, FlowRate, FlowRateLike};
use crate::utils::{all_same_value, unwrap_all};
use itertools::Itertools;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

//
// Edge struct for Compacted Edbg
//

/// Basic implementations of Compacted EDbgEdge
#[derive(Clone, Debug)]
pub struct SimpleCompactedEDbgEdge<K: KmerLike> {
    kmer: K,
    copy_num: CopyNum,
    origin_edges: Vec<EdgeIndex>,
    // origin_nodes: Vec<NodeIndex>,
}

impl<K: KmerLike> SimpleCompactedEDbgEdge<K> {
    pub fn new(kmer: K, copy_num: CopyNum, origin_edges: Vec<EdgeIndex>) -> Self {
        SimpleCompactedEDbgEdge {
            kmer,
            copy_num,
            origin_edges,
        }
    }
    pub fn origin_edges(&self) -> &[EdgeIndex] {
        &self.origin_edges
    }
    pub fn copy_num(&self) -> CopyNum {
        self.copy_num
    }
}

impl<K: KmerLike> std::fmt::Display for SimpleCompactedEDbgEdge<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{} (x{}) ({})",
            self.kmer,
            self.copy_num,
            self.origin_edges.iter().map(|e| e.index()).join(",")
        )
    }
}

///
/// Flow in Compacted-EDbg -> Flow in EDbg (Flow is EdgeVec<DenseStorage>)
/// single edge in compacted-edbg corresponds to multiple edges in edbg.
///
/// * n_kmers: total number of kmers in original edbg
///     = n_nodes in Dbg
///     = n_edges in (uncompacted) EDbg
///
pub fn compacted_flow_into_original_flow<N, K, F>(
    n_kmers: usize,
    compacted_edbg: &DiGraph<N, SimpleCompactedEDbgEdge<K>>,
    compacted_flow: &Flow<F>,
) -> Flow<F>
where
    K: KmerLike,
    F: FlowRateLike,
{
    let mut ret: Vec<Option<F>> = vec![None; n_kmers];
    for edge in compacted_edbg.edge_indices() {
        let weight = compacted_edbg.edge_weight(edge).unwrap();
        let f = compacted_flow[edge];
        for &origin_edge in weight.origin_edges() {
            ret[origin_edge.index()] = Some(f);
        }
    }
    Flow::from_inner_vec(unwrap_all(ret))
}

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    /// create simple-path collapsed edbg
    ///
    pub fn to_compact_edbg_graph(
        &self,
    ) -> DiGraph<SimpleEDbgNode<N::Kmer>, SimpleCompactedEDbgEdge<N::Kmer>> {
        let graph: DiGraph<_, _> = self.to_edbg_graph(
            |km1mer| SimpleEDbgNode::new(km1mer.clone()),
            |node, weight| {
                let kmer = weight.kmer().clone();
                let copy_num = weight.copy_num();
                SimpleEDbgEdge::new(kmer, copy_num, node)
            },
        );
        // println!("{}", petgraph::dot::Dot::with_config(&graph, &[]));
        let compacted = compact_simple_paths(&graph);
        compacted.map(
            |_node, weight| weight.clone(),
            |_edge, weight| {
                let edges: Vec<_> = weight.iter().map(|(e, _)| *e).collect();
                let kmer = weight
                    .iter()
                    .map(|(_, w)| w.kmer().clone())
                    .reduce(|accum, kmer| accum.overlap(&kmer))
                    .unwrap();
                let copy_num = all_same_value(weight.iter().map(|(_, w)| w.copy_num()))
                    .expect("not all copynums in edge are the same");
                SimpleCompactedEDbgEdge::new(kmer, copy_num, edges)
            },
        )
    }
}
