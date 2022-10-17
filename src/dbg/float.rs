//!
//! de bruijn graph with float (real-valued) copy numbers
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase};
use crate::common::Seq;
use crate::graph::float_seq_graph::{FloatSeqEdge, FloatSeqGraph, FloatSeqNode};
use crate::hmmv2::common::PModel;
use crate::hmmv2::q::QScore;
use crate::hmmv2::{EdgeFreqs, NodeFreqs};
use crate::min_flow::FlowRateLike;
use crate::prelude::*;
use crate::vector::{DenseStorage, EdgeVec, NodeVec};
use fnv::FnvHashMap as HashMap;
use itertools::{iproduct, izip, Itertools};
use petgraph::dot::Dot;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use rayon::prelude::IntoParallelIterator;

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
    pub fn scale_by_total_density(&mut self, total_density: CopyDensity) {
        self.scale_density(total_density / self.total_density());
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
    /// get the number of redundant nodes, that is the copy density is below the specified min_density.
    ///
    pub fn n_redundant_nodes(&self, min_density: CopyDensity) -> usize {
        self.nodes()
            .filter(|(_, weight)| weight.copy_density() < min_density)
            .count()
    }
    ///
    /// get the number of dead (= x0) nodes.
    ///
    pub fn n_dead_nodes(&self) -> usize {
        self.nodes()
            .filter(|(_, weight)| weight.copy_density() == 0.0)
            .count()
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
        let (_, parent, _) = self.parents(node).next().unwrap_or_else(|| {
            panic!(
                "node {}(#{}) has no parent (copy_density={})",
                self.node(node).kmer(),
                node.index(),
                self.node(node).copy_density(),
            )
        });
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
    ///
    pub fn to_full_prob_parallel<T>(&self, param: PHMMParams, seqs: T) -> Prob
    where
        T: IntoParallelIterator,
        T::Item: Seq,
    {
        self.to_phmm(param).to_full_prob_parallel(seqs)
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
    ///
    ///
    pub fn set_node_copy_densities(&mut self, copy_densities: &NodeVec<DenseStorage<CopyDensity>>) {
        assert!(copy_densities.len() == self.n_nodes());
        for (i, node_weight_mut) in self.graph.node_weights_mut().enumerate() {
            let node = NodeIndex::new(i);
            let copy_density = copy_densities[node];
            node_weight_mut.set_copy_density(copy_density);
        }
    }
}

//
// dbg consistency validations
//
impl<K: KmerLike> Dbg<FloatDbgNode<K>, FloatDbgEdge> {
    pub fn is_valid(&self) -> bool {
        self.is_graph_valid() && self.has_no_duplicated_node() && self.has_no_parallel_edge()
    }
}

//
// k->k+1 upgrade related
//
impl<K: KmerLike> Dbg<FloatDbgNode<K>, FloatDbgEdge> {
    ///
    /// get a naive copy density of an edge
    ///
    /// ```text
    /// copy_density(v->w) = copy_density(v) * (copy_density(w) / copy_density(v's childs))
    /// ```
    ///
    fn naive_edge_copy_density(&self, edge: EdgeIndex) -> CopyDensity {
        let (v, w) = self.edge_endpoints(edge).expect("edge does not exist");

        // sum of copy_density of childs of v
        let dws: CopyDensity = self
            .childs(v)
            .map(|(_, w, _)| self.node(w).copy_density())
            .sum();
        let dv = self.node(v).copy_density();
        let dw = self.node(w).copy_density();

        if dws == 0.0 {
            0.0
        } else {
            dv * (dw / dws)
        }
    }
    ///
    /// Create a `k+1` dbg from the `k` dbg using `naive_edge_copy_density`
    ///
    pub fn to_kp1_dbg(&self) -> Self {
        let mut graph = DiGraph::new();

        // mapping from "edge in k-dbg" into "node in k+1-dbg".
        let mut ids: HashMap<EdgeIndex, NodeIndex> = HashMap::default();

        // (1) nodes/edges outside tip area
        // a edge in k-dbg is corresponds to a node in k+1-dbg.
        for (edge, s, t, _weight) in self.edges() {
            if !self.is_warp_edge(edge) {
                let kmer = self.kmer(s).join(self.kmer(t));
                let copy_density = self.naive_edge_copy_density(edge);
                let node = graph.add_node(FloatDbgNode::new(kmer, copy_density));
                ids.insert(edge, node);
            }
        }
        for (node, weight) in self.nodes() {
            if !weight.is_head() && !weight.is_tail() {
                // add an edge between all in-edges and out-edges pair.
                // --e1--> node --e2--> in k-dbg
                for (e1, _, _) in self.parents(node) {
                    let v1 = ids.get(&e1).unwrap();
                    for (e2, _, _) in self.childs(node) {
                        let v2 = ids.get(&e2).unwrap();
                        // copy numbers of edges in k+1 is ambiguous.
                        graph.add_edge(*v1, *v2, FloatDbgEdge::new(None));
                    }
                }
            }
        }

        // (2) nodes/edges in tip area
        // intersections
        let tips = self.tips_base();
        let in_nodes: Vec<NodeIndex> = tips
            .iter_in_node_indexes()
            .map(|v| {
                // add a node of in_node
                graph.add_node(FloatDbgNode::new(
                    self.node(v).kmer().extend_tail(),
                    self.node(v).copy_density(),
                ))
            })
            .collect();
        let out_nodes: Vec<NodeIndex> = tips
            .iter_out_node_indexes()
            .map(|v| {
                // add a node of out_node
                graph.add_node(FloatDbgNode::new(
                    self.node(v).kmer().extend_head(),
                    self.node(v).copy_density(),
                ))
            })
            .collect();
        for (&w1, &w2) in iproduct!(in_nodes.iter(), out_nodes.iter()) {
            graph.add_edge(w1, w2, FloatDbgEdge::new(None));
        }

        // add an edge for a in_edges of tail nodes
        for (v, &w) in izip!(tips.iter_in_node_indexes(), in_nodes.iter()) {
            for (e, _, _) in self.parents(v) {
                let w2 = ids.get(&e).unwrap();
                // w2: parent(YXNN) -> w: tail(XNNN)
                graph.add_edge(*w2, w, FloatDbgEdge::new(None));
            }
        }

        // add an edge for a out_edges of head nodes
        for (v, &w) in izip!(tips.iter_out_node_indexes(), out_nodes.iter()) {
            for (e, _, _) in self.childs(v) {
                let w2 = ids.get(&e).unwrap();
                // w: head(NNNX) -> w2: child(NNXY)
                graph.add_edge(w, *w2, FloatDbgEdge::new(None));
            }
        }

        Self::from_digraph(self.k() + 1, graph)
    }
    ///
    /// remove nodes whose copy density is 0.0
    ///
    pub fn remove_zero_copy_node(&mut self) {
        // println!("removing {} nodes", self.n_dead_nodes());
        // TODO rounding error causes small nodes that have infinitesimal copy density. this
        // function removes these nodes.
        self.graph
            .retain_nodes(|g, v| g.node_weight(v).unwrap().copy_density() >= f64::eps());
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
    fn float_dbg_kp1_convert() {
        let mut dbg = mock_intersection();
        let mut fdbg = FloatDbg::from_dbg(&dbg);
        println!("{}", fdbg);
        let fdbg_kp1 = fdbg.to_kp1_dbg();
        println!("{}", fdbg_kp1);
        // for (k, c) in sorted_node_list(&fdbg_kp1) {
        //     println!("{} {}", k, c);
        // }
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
