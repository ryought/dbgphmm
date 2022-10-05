//!
//! de bruijn graph with float (real-valued) copy numbers
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase};
use crate::graph::float_seq_graph::{FloatSeqEdge, FloatSeqGraph, FloatSeqNode};
use crate::hmmv2::common::PModel;
use crate::hmmv2::q::QScore;
use crate::hmmv2::{EdgeFreqs, NodeFreqs};
use crate::prelude::*;
use crate::vector::{DenseStorage, EdgeVec, NodeVec};
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
    pub fn new(kmer: K, copy_density: CopyDensity) -> Self {
        FloatDbgNode { kmer, copy_density }
    }
    pub fn copy_density(&self) -> CopyDensity {
        self.copy_density
    }
    pub fn round_copy_num(&self) -> CopyNum {
        self.copy_density.round() as CopyNum
    }
    pub fn set_copy_density(&mut self, copy_density: CopyDensity) {
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
// Copy density and PHMM related
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
    /// modify node's copy_density by calling `FloatDbgNode::set_copy_density` without checking
    /// flow consistency.
    ///
    pub fn set_node_copy_density(&mut self, node: NodeIndex, copy_density: CopyDensity) {
        self.graph
            .node_weight_mut(node)
            .unwrap()
            .set_copy_density(copy_density)
    }
    ///
    /// calculate the total density of all nodes
    ///
    pub fn total_density(&self) -> CopyDensity {
        self.nodes().map(|(v, vw)| vw.copy_density()).sum()
    }
    ///
    /// calculate the total density of all emittable nodes
    ///
    pub fn total_emittable_copy_density(&self) -> CopyDensity {
        self.graph.total_emittable_copy_density()
    }
    ///
    /// calculate the density of nodes in the intersection (with the given node)
    ///
    pub fn intersection_emittable_copy_density(&self, node: NodeIndex) -> CopyDensity {
        let (_, parent, _) = self.parents(node).next().unwrap();
        self.childs(parent)
            .filter(|(_, sibling, _)| self.node(*sibling).is_emittable())
            .map(|(_, sibling, _)| self.node(sibling).copy_density)
            .sum()
    }
    ///
    /// convert to PHMM `PModel` by running `FloatSeqGraph::to_phmm` for `self.graph`.
    ///
    pub fn to_phmm(&self, param: PHMMParams) -> PModel {
        self.graph.to_phmm(param)
    }
    ///
    ///
    pub fn to_node_copy_densities(&self) -> NodeVec<DenseStorage<CopyDensity>> {
        let mut v = NodeVec::new(self.n_nodes(), 0.0);
        for (node, weight) in self.nodes() {
            v[node] = weight.copy_density();
        }
        v
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
/// q_score_diff_exact = init + trans
/// init = A + B
/// trans = C + D
/// A = A[0,l] * (log(c[l]+diff) - log(c[l]))
/// B = - sum(A[0,l]) * (log(G+diff) - log(G))
/// C = sum(A[k,l]) * (log(c[l]+diff) - log(c[l]))
/// D = - sum(A[i]) * (log(G[i]+diff) - log(G[i]))
/// ```
///
pub fn q_score_diff_exact<K: KmerLike>(
    dbg: &FloatDbg<K>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    node: NodeIndex,
    diff: CopyDensity,
) -> QScore {
    let copy_density = dbg.node(node).copy_density;
    let copy_density_total: CopyDensity = dbg.total_emittable_copy_density();
    let copy_density_intersection: CopyDensity = dbg.intersection_emittable_copy_density(node);
    let is_emittable = dbg.node(node).is_emittable();

    // init
    // A: A[0,l]
    let init = if is_emittable {
        let a0: CopyDensity = init_freqs[node];
        let a = a0 * ln_diff(copy_density, diff);
        // B: G
        let b0: CopyDensity = dbg
            .nodes()
            .filter(|(_, vw)| vw.is_emittable())
            .map(|(v, _)| init_freqs[v])
            .sum();
        let b = -b0 * ln_diff(copy_density_total, diff);
        a + b
    } else {
        0.0
    };

    // trans
    // C: sum(A[k,l])
    let trans = if is_emittable {
        let c0: CopyDensity = dbg
            .parents(node)
            .filter(|(_, parent, _)| dbg.is_emittable(*parent))
            .map(|(e, _, _)| edge_freqs[e])
            .sum();
        let c = c0 * ln_diff(copy_density, diff);
        // D: Gi
        let d0: CopyDensity = dbg
            .parents(node)
            .filter(|(_, parent, _)| dbg.is_emittable(*parent))
            .map(|(_, parent, _)| {
                dbg.childs(parent)
                    .filter(|(_, child, _)| dbg.is_emittable(*child))
                    .map(|(e, child, _)| edge_freqs[e])
                    .sum::<CopyDensity>()
            })
            .sum();
        let d = -d0 * ln_diff(copy_density_intersection, diff);
        c + d
    } else {
        0.0
    };

    QScore::new(init, trans, 0.0)
}

///
/// log_e(x + dx) - log_e(x)
///
/// TODO is there any better/precise/stable way to calculate the value?
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
    use crate::hmmv2::params::PHMMParams;
    use crate::hmmv2::q::q_score_exact;

    #[test]
    fn convert_to_float_dbg() {
        let mut dbg = mock_intersection_small();
        println!("{}", dbg);
        let mut fdbg = FloatDbg::from_dbg(&dbg);
        fdbg.scale_density(0.2);
        println!("{}", fdbg);
        println!("td={}", fdbg.total_density());
    }

    #[test]
    fn compare_q_score_and_q_score_diff() {
        //
        let mut dbg = mock_intersection_small();
        let mut fdbg = FloatDbg::from_dbg(&dbg);
        println!("{}", fdbg);
        let reads = [b"ATAGCT"];

        // normal score
        let phmm = fdbg.to_phmm(PHMMParams::zero_error());
        let (edge_freqs, init_freqs, full_prob) = phmm.to_edge_and_init_freqs_parallel(&reads);
        let q0 = q_score_exact(&phmm, &edge_freqs, &init_freqs);
        println!("{}", q0);

        // density with multiplying constant should result in the same score
        fdbg.scale_density(100.0);
        let phmm = fdbg.to_phmm(PHMMParams::zero_error());
        let (edge_freqs, init_freqs, full_prob) = phmm.to_edge_and_init_freqs_parallel(&reads);
        let q = q_score_exact(&phmm, &edge_freqs, &init_freqs);
        assert_abs_diff_eq!(q0.total(), q.total(), epsilon = 0.00001);
        println!("{}", q);

        // density with multiplying constant should result in the same score
        let diff = 10.0;
        for (node, _) in dbg.nodes() {
            println!("{:?}", node);
            // (1) modify phmm
            let mut fdbg_mod = fdbg.clone();
            fdbg_mod.set_node_copy_density(node, fdbg.node(node).copy_density + diff);
            let phmm = fdbg_mod.to_phmm(PHMMParams::zero_error());
            let q1 = q_score_exact(&phmm, &edge_freqs, &init_freqs);
            let q1d = q1.sub(q);
            // println!("{} {}", fdbg_mod, phmm);

            // (2) calculate only diff
            let q2d = q_score_diff_exact(&fdbg, &edge_freqs, &init_freqs, node, diff);
            println!("{} {}", q1d, q2d);
            assert_abs_diff_eq!(q1d.total(), q2d.total());
        }
    }
}
