//!
//! Float dbg optimization by EM & GradDescent
//!
//! * takes FloatDbg as input
//! * E-step: run forward/backward algorithm and calculate the node/edge usage
//! * M-step: improve Q score by GradDescent and MinCostFlow
//! * iterate E/M-steps and returns the improved FloatDbg
//!
use crate::dbg::dbg::DbgNodeBase;
use crate::dbg::edge_centric::EDbgEdgeBase;
use crate::dbg::float::{q_score_diff_exact, CopyDensity, FloatDbg, FloatDbgEdge, FloatDbgNode};
use crate::graph::float_seq_graph::FloatSeqGraph;
use crate::hmmv2::{EdgeFreqs, NodeFreqs};
use crate::min_flow::residue::{
    improve_residue_graph, ResidueDirection, ResidueEdge, ResidueGraph,
};
use crate::prelude::*;

///
/// E-step: calculate edge_freqs (freq between v->w) and init_freqs (freq between Begin->w)
///
pub fn e_step<K: KmerLike>(
    dbg: &FloatDbg<K>,
    reads: &Reads,
    params: &PHMMParams,
) -> (EdgeFreqs, NodeFreqs, Prob) {
    let phmm = dbg.graph.to_phmm(params.clone());
    phmm.to_edge_and_init_freqs_parallel(reads)
}

///
/// M-step
///
pub fn m_step<K: KmerLike>(
    dbg: &FloatDbg<K>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyDensity,
) {
    // (1) convert to edge-centric dbg with each edge has a cost

    // (2) search for negative cycle

    // search again for
}

///
/// FloatDbg
/// -> construct residue graph of QScore difference when +1/-1
/// -> change density if negative meaningful cycle was found
/// -> new updated FloatDbg
///
pub fn m_step_once<K: KmerLike>(
    dbg: &FloatDbg<K>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyDensity,
) -> FloatDbg<K> {
    unimplemented!();
}

#[derive(Clone, Debug)]
struct QDiffEdge<K: KmerLike> {
    // q_inc: f64,
    // q_dec: f64,
    origin_node: NodeIndex,
    copy_density: CopyDensity,
    kmer: K,
}
impl<K: KmerLike> QDiffEdge<K> {
    fn new(origin_node: NodeIndex, copy_density: CopyDensity, kmer: K) -> Self {
        QDiffEdge {
            origin_node,
            copy_density,
            kmer,
        }
    }
}
impl<K: KmerLike> EDbgEdgeBase for QDiffEdge<K> {
    type Kmer = K;
    fn kmer(&self) -> &K {
        &self.kmer
    }
    fn origin_node(&self) -> NodeIndex {
        self.origin_node
    }
}

///
/// convert to residue graph, each edge has q_score_diff of +/- diff.
///
/// we want to maximize QScore and MinCostFlow solves minimize the cost, so the negated q_score_diff is used as a edge cost.
///
fn to_residue_graph<K: KmerLike>(
    dbg: &FloatDbg<K>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    diff: CopyDensity,
) -> ResidueGraph {
    let edbg = dbg.to_edbg_generic(
        |_kmer| (),
        |node, node_weight| {
            QDiffEdge::new(node, node_weight.copy_density(), node_weight.kmer().clone())
        },
    );
    let max_copy_density = 10000.0;
    let mut rg: ResidueGraph = ResidueGraph::new();
    for (e, v, w, ew) in edbg.edges() {
        let node = ew.origin_node();
        let mut edges = Vec::new();
        if ew.copy_density < max_copy_density {
            // increasable
            edges.push((
                v,
                w,
                ResidueEdge::new(
                    1,
                    -q_score_diff_exact(dbg, edge_freqs, init_freqs, node, diff).total(),
                    e,
                    ResidueDirection::Up,
                ),
            ));
        }
        if ew.copy_density > diff {
            // decreasable
            edges.push((
                w,
                v,
                ResidueEdge::new(
                    1,
                    -q_score_diff_exact(dbg, edge_freqs, init_freqs, node, -diff).total(),
                    e,
                    ResidueDirection::Down,
                ),
            ));
        }
        rg.extend_with_edges(&edges);
    }
    rg
}

fn apply_to_dbg<K: KmerLike>(
    dbg: &mut FloatDbg<K>,
    diff: CopyDensity,
    rg: &ResidueGraph,
    edges: &[EdgeIndex],
) {
    for e in edges {
        let ew = rg.edge_weight(*e).unwrap();
        let node = NodeIndex::new(ew.target.index());
        let copy_density = dbg.node(node).copy_density();
        let new_copy_density = match ew.direction {
            ResidueDirection::Up => copy_density + diff,
            ResidueDirection::Down => copy_density - diff,
        };
        dbg.set_node_copy_density(node, new_copy_density);
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks::*;
    use petgraph::dot::Dot;

    #[test]
    fn em_float_residue_graph() {
        let dbg = mock_intersection_small();
        let reads = [b"ATAGCT"];
        let mut fdbg = FloatDbg::from_dbg(&dbg);
        let phmm = fdbg.to_phmm(PHMMParams::zero_error());
        let (edge_freqs, init_freqs, full_prob) = phmm.to_edge_and_init_freqs_parallel(&reads);
        let diff = 0.1;
        let rg = to_residue_graph(&fdbg, &edge_freqs, &init_freqs, diff);
        println!("{:?}", Dot::with_config(&rg, &[]));
        let edges = improve_residue_graph(&rg).unwrap();
        apply_to_dbg(&mut fdbg, diff, &rg, &edges);
        println!("{}", fdbg);
    }
}
