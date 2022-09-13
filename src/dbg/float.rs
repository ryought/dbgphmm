//!
//! de bruijn graph with float (real-valued) copy numbers
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase};
use crate::graph::float_seq_graph::{FloatSeqEdge, FloatSeqNode};
use crate::hmmv2::q::QScore;
use crate::hmmv2::{EdgeFreqs, NodeFreqs};
use crate::prelude::*;
use petgraph::dot::Dot;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

/// `CopyDensity` = f64
/// Float valued copy number
pub type CopyDensity = f64;

/// de bruijn graph with float (real-valued) copy numbers
pub type FloatDbg<K> = Dbg<FloatDbgNode<K>, FloatDbgEdge>;

/// node struct of FloatDbg
#[derive(Clone)]
pub struct FloatDbgNode<K: KmerLike> {
    kmer: K,
    copy_density: CopyDensity,
}

impl<K: KmerLike> DbgNodeBase for FloatDbgNode<K> {
    type Kmer = K;
    fn kmer(&self) -> &K {
        &self.kmer
    }
}

impl<K: KmerLike> FloatDbgNode<K> {
    fn new(kmer: K, copy_density: CopyDensity) -> Self {
        FloatDbgNode { kmer, copy_density }
    }
    fn copy_density(&self) -> CopyDensity {
        self.copy_density
    }
    fn round_copy_num(&self) -> CopyNum {
        self.copy_density.round() as CopyNum
    }
    fn set_copy_density(&mut self, copy_density: CopyDensity) {
        self.copy_density = copy_density
    }
}

/// edge struct of FloatDbg
#[derive(Clone)]
pub struct FloatDbgEdge {
    copy_density: Option<CopyDensity>,
}

impl FloatDbgEdge {
    fn new(copy_density: Option<CopyDensity>) -> Self {
        FloatDbgEdge { copy_density }
    }
    fn copy_density(&self) -> Option<CopyDensity> {
        self.copy_density
    }
    fn round_copy_num(&self) -> Option<CopyNum> {
        match self.copy_density {
            Some(copy_density) => Some(copy_density.round() as CopyNum),
            None => None,
        }
    }
    fn set_copy_density(&mut self, copy_density: Option<CopyDensity>) {
        self.copy_density = copy_density
    }
}

//
// Constructors
//
impl<K: KmerLike> Dbg<FloatDbgNode<K>, FloatDbgEdge> {
    ///
    /// create FloatDbg from (normal, integer copy numbered) Dbg.
    ///
    pub fn from_dbg<N: DbgNode<Kmer = K>, E: DbgEdge>(dbg: &Dbg<N, E>) -> Self {
        let g = dbg.graph.map(
            |_, vw| FloatDbgNode::new(vw.kmer().clone(), vw.copy_num() as CopyDensity),
            |_, ew| FloatDbgEdge::new(ew.copy_num().map(|copy_num| copy_num as CopyDensity)),
        );
        FloatDbg::from_digraph(dbg.k(), g)
    }
}

//
// Copy density related
//
impl<K: KmerLike> Dbg<FloatDbgNode<K>, FloatDbgEdge> {
    ///
    /// scale all copy densities by constant factor
    /// i.e. multiply the given `scale` to all node/edge copy_density(s).
    ///
    pub fn scale_density(&mut self, scale: CopyDensity) {
        for nw in self.graph.node_weights_mut() {
            let new_copy_density = nw.copy_density() * scale;
            nw.set_copy_density(new_copy_density);
        }
        for ew in self.graph.edge_weights_mut() {
            let new_copy_density = ew.copy_density().map(|d| d * scale);
            ew.set_copy_density(new_copy_density);
        }
    }
    ///
    /// calculate the total density of all nodes
    ///
    pub fn total_density(&self) -> CopyDensity {
        self.nodes().map(|(v, vw)| vw.copy_density()).sum()
    }
}

//
// std::fmt::Display
//
impl<K: KmerLike> std::fmt::Display for FloatDbgNode<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} (x{:.5})", self.kmer(), self.copy_density())
    }
}
impl std::fmt::Display for FloatDbgEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.copy_density() {
            Some(copy_density) => write!(f, "x{:.5}", copy_density),
            None => write!(f, "x?"),
        }
    }
}
impl<K: KmerLike> std::fmt::Display for FloatDbg<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.graph, &[]))
    }
}

//
// FloatSeqGraph
// to convert PHMMModel
//
impl<K: KmerLike> FloatSeqNode for FloatDbgNode<K> {
    fn copy_density(&self) -> CopyDensity {
        self.copy_density()
    }
    fn base(&self) -> u8 {
        self.emission()
    }
}
impl FloatSeqEdge for FloatDbgEdge {
    fn copy_density(&self) -> Option<CopyDensity> {
        self.copy_density()
    }
}

//
// Q scores
//
///
/// Calculate difference of QScore
/// when `node`'s copy density was changed by `diff`.
///
/// ```text
/// q_score_diff_exact = A + B
/// A = (sum(A[0,l]) + sum(A[k,l])) * (log(c[l]+diff) - log(c[l]))
/// B = - sum(A[i]) * (log(G[i]+diff) - log(G[i]))
/// ```
///
pub fn q_score_diff_exact<K: KmerLike>(
    dbg: &FloatDbg<K>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    node: NodeIndex,
    diff: CopyDensity,
) -> QScore {
    // A
    let a0 = init_freqs[node];
    let a1 = edge_freqs;
    unimplemented!();
}

///
/// log_e(x + dx) - log_e(x)
///
pub fn ln_diff(x: f64, dx: f64) -> f64 {
    (x + dx).ln() - x.ln()
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::super::mocks::*;
    use super::*;

    #[test]
    fn convert_to_float_dbg() {
        let mut dbg = mock_intersection_small();
        println!("{}", dbg);
        let mut fdbg = FloatDbg::from_dbg(&dbg);
        fdbg.scale_density(0.2);
        println!("{}", fdbg);
        println!("td={}", fdbg.total_density());
    }
}
