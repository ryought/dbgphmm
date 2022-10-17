//!
//! Float dbg optimization by EM & GradDescent
//!
//! * takes FloatDbg as input
//! * E-step: run forward/backward algorithm and calculate the node/edge usage
//! * M-step: improve Q score by GradDescent and MinCostFlow
//! * iterate E/M-steps and returns the improved FloatDbg
//!
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase};
use crate::dbg::edge_centric::EDbgEdgeBase;
use crate::dbg::float::{q_score_diff_exact, CopyDensity, FloatDbg, FloatDbgEdge, FloatDbgNode};
use crate::graph::float_seq_graph::FloatSeqGraph;
use crate::hist::stat;
use crate::hmmv2::q::{q_score_exact, QScore};
use crate::hmmv2::{EdgeFreqs, NodeFreqs};
use crate::io::cytoscape::{NodeAttr, NodeAttrVec};
use crate::min_flow::flow::{ConstCost, Flow, FlowEdge};
use crate::min_flow::min_cost_flow_from;
use crate::min_flow::residue::{
    improve_residue_graph, ResidueDirection, ResidueEdge, ResidueGraph,
};
use crate::prelude::*;
use crate::vector::{DenseStorage, NodeVec};

///
/// run e-step and m-step iteratively
///
pub fn em<K: KmerLike>(
    dbg: &FloatDbg<K>,
    reads: &Reads,
    // parameters
    params: &PHMMParams,
    genome_size: CopyDensity,
    diff: CopyDensity,
    n_max_em_iteration: usize,
    n_max_iteration: usize,
) -> EMResult<K> {
    let mut dbg_current = dbg.clone();
    let mut ret = EMResult::new();

    for i in 0..n_max_em_iteration {
        let (edge_freqs, init_freqs, p) = e_step(&dbg_current, &reads, &params);
        eprintln!("#{} p={}", i, p);
        ret.e.push(p);

        let (dbg_new, m_step_result) = m_step(
            &dbg_current,
            &edge_freqs,
            &init_freqs,
            genome_size,
            diff,
            n_max_iteration,
            &params,
        );
        ret.m.push(m_step_result);

        match dbg_new {
            Some(dbg_new) => {
                dbg_current = dbg_new;
            }
            None => break,
        };
    }

    ret
}

pub struct EMResult<K: KmerLike> {
    pub e: Vec<Prob>,
    pub m: Vec<Vec<MStepResult<K>>>,
}

impl<K: KmerLike> EMResult<K> {
    pub fn new() -> Self {
        EMResult {
            e: Vec::new(),
            m: Vec::new(),
        }
    }
}

///
///
pub fn em_result_to_final_dbg<K: KmerLike>(result: &EMResult<K>) -> Option<FloatDbg<K>> {
    for (em_id, m_step_result) in result.m.iter().rev().enumerate() {
        for (m_id, m_step_once_result) in m_step_result.iter().rev().enumerate() {
            match m_step_once_result {
                MStepResult::Update(dbg, _) => {
                    return Some(dbg.clone());
                }
                _ => {}
            };
        }
    }
    return None;
}

///
/// create Vec<NodeAttrVec> (that represents the time series of copy densities of nodes of EM
/// steps) by using true dbg (Dbg<N, E>, not floated) and result (EMResult<K>).
///
pub fn em_result_to_node_historys<K: KmerLike>(
    result: &EMResult<K>,
) -> Vec<(String, NodeVec<DenseStorage<CopyDensity>>)> {
    let mut ret = Vec::new();

    for (em_id, m_step_result) in result.m.iter().enumerate() {
        for (m_id, m_step_once_result) in m_step_result.iter().enumerate() {
            match m_step_once_result {
                MStepResult::Update(dbg, _) => {
                    let copy_densities = dbg.to_node_copy_densities();
                    let label = format!("{}#{}", em_id, m_id);
                    ret.push((label, copy_densities));
                }
                _ => {}
            }
        }
    }

    ret
}

///
/// inspect the relationship between density and copy number
///
/// for each copy number c, draw the histogram of density of nodes whose copy number is c.
///
pub fn inspect_density_histgram<N: DbgNode, E: DbgEdge, K: KmerLike>(
    dbg_true: &Dbg<N, E>,
    final_fdbg: &FloatDbg<K>,
) {
    let densitys = final_fdbg.to_node_copy_densities();
    inspect_freqs_histgram(dbg_true, &densitys)
}

///
/// inspect the relationship between density and copy number
///
/// for each copy number c, draw the histogram of density of nodes whose copy number is c.
///
pub fn inspect_freqs_histgram<N: DbgNode, E: DbgEdge>(
    dbg_true: &Dbg<N, E>,
    densitys: &NodeVec<DenseStorage<CopyDensity>>,
) {
    let copy_nums_list = dbg_true.to_copy_nums_list();
    for (copy_num, nodes) in copy_nums_list {
        let densitys_of_copy_num: Vec<_> = nodes.iter().map(|&node| densitys[node]).collect();
        let (ave, std, min, max) = stat(&densitys_of_copy_num);
        eprintln!(
            "[{}x] ave={} std={} min={} max={}",
            copy_num, ave, std, min, max
        );
        // eprintln!("{:?}", densitys);
    }
}

enum ShrinkEdgeType {
    Redundant,
    NotRedundant,
}

struct ShrinkEdge {
    current_density: CopyDensity,
    edge_type: ShrinkEdgeType,
}

impl FlowEdge<f64> for ShrinkEdge {
    fn demand(&self) -> f64 {
        match self.edge_type {
            ShrinkEdgeType::Redundant => 0.0,
            ShrinkEdgeType::NotRedundant => self.current_density,
        }
    }
    fn capacity(&self) -> f64 {
        match self.edge_type {
            ShrinkEdgeType::Redundant => self.current_density,
            ShrinkEdgeType::NotRedundant => 1000.0, // TODO some large constant
        }
    }
}

impl ConstCost for ShrinkEdge {
    fn cost(&self) -> f64 {
        match self.edge_type {
            ShrinkEdgeType::Redundant => 1.0,
            ShrinkEdgeType::NotRedundant => 0.0,
        }
    }
}

///
/// shrink nodes whose `copy_density` is less than `min_density` with satisfying flow constaint.
///
pub fn shrink_nodes<K: KmerLike>(
    fdbg: &FloatDbg<K>,
    min_density: CopyDensity,
) -> NodeVec<DenseStorage<CopyDensity>> {
    // convert to flow network
    let edbg = fdbg.to_edbg_generic(
        // to_node:
        |_| (),
        // to_edge:
        |node, weight| {
            let density = weight.copy_density();
            ShrinkEdge {
                current_density: density,
                edge_type: if density < min_density {
                    ShrinkEdgeType::Redundant
                } else {
                    ShrinkEdgeType::NotRedundant
                },
            }
        },
    );
    // solve as min_flow
    let current_flow: Flow<f64> = fdbg.to_node_copy_densities().switch_index();
    let flow = min_cost_flow_from(&edbg.graph, &current_flow);
    flow.switch_index()
}

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
    diff: CopyDensity,
    n_max_iteration: usize,
    params: &PHMMParams,
) -> (Option<FloatDbg<K>>, Vec<MStepResult<K>>) {
    let mut dbg_current = dbg.clone();
    let mut is_updated = false;
    let mut ret = Vec::new();

    for i in 0..n_max_iteration {
        let r = m_step_once(
            &dbg_current,
            edge_freqs,
            init_freqs,
            genome_size,
            diff,
            params,
        );

        ret.push(r.clone());

        match r {
            MStepResult::Update(dbg_new, q_score) => {
                eprintln!("update {}", q_score);
                dbg_current = dbg_new;
                is_updated = true;
            }
            MStepResult::NoImprove(dbg_new, q_score) => {
                eprintln!("no improve {}", q_score);
                break;
            }
            MStepResult::NoNegCycle => {
                eprintln!("no neg cycle");
                break;
            }
        };
    }

    if is_updated {
        // some new dbg was generated by this m_step
        (Some(dbg_current), ret)
    } else {
        // this m_step did not change dbg
        (None, ret)
    }
}

#[derive(Clone)]
pub enum MStepResult<K: KmerLike> {
    /// The found negative cycle improved q-score
    Update(FloatDbg<K>, QScore),
    /// A negative cycle was found, but it did not improve q-score
    NoImprove(FloatDbg<K>, QScore),
    /// Any negative cycle was not found.
    NoNegCycle,
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
    diff: CopyDensity,
    params: &PHMMParams,
) -> MStepResult<K> {
    //(0) calculate the original qscore
    let q_score = q_score_exact(&dbg.to_phmm(params.clone()), edge_freqs, init_freqs);

    let mut fdbg = dbg.clone();
    // (1) convert to edge-centric dbg with each edge has a cost
    let rg = to_residue_graph(&fdbg, &edge_freqs, &init_freqs, diff);
    // (2) search for negative cycle
    match improve_residue_graph(&rg) {
        Some(edges) => {
            apply_to_dbg(&mut fdbg, diff, &rg, &edges);
            let q_score_new = q_score_exact(&fdbg.to_phmm(params.clone()), edge_freqs, init_freqs);

            // (3) check if the copy density changes actually improves q-score.
            if q_score_new.total() >= q_score.total() {
                MStepResult::Update(fdbg, q_score_new)
            } else {
                MStepResult::NoImprove(fdbg, q_score_new)
            }
        }
        None => MStepResult::NoNegCycle,
    }
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
) -> ResidueGraph<usize> {
    let edbg = dbg.to_edbg_generic(
        |_kmer| (),
        |node, node_weight| {
            QDiffEdge::new(node, node_weight.copy_density(), node_weight.kmer().clone())
        },
    );
    let max_copy_density = 10000.0;
    let mut rg = ResidueGraph::new();
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
    rg: &ResidueGraph<usize>,
    edges: &[EdgeIndex],
) {
    let genome_size_orig = dbg.total_density();

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

    // scale copy densities of all nodes
    // so that this function does not change the genome size
    let genome_size_new = dbg.total_density();
    let scale = genome_size_orig / genome_size_new;
    dbg.scale_density(scale);
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks::*;
    use crate::kmer::veckmer::kmer;
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

        assert_abs_diff_eq!(
            fdbg.node(fdbg.find_node_from_kmer(&kmer(b"AGCT")).unwrap())
                .copy_density(),
            1.1,
        );
        assert_abs_diff_eq!(
            fdbg.node(fdbg.find_node_from_kmer(&kmer(b"AGCC")).unwrap())
                .copy_density(),
            0.9,
        );
        assert_abs_diff_eq!(
            fdbg.node(fdbg.find_node_from_kmer(&kmer(b"nnnA")).unwrap())
                .copy_density(),
            1.0,
        );
        assert_abs_diff_eq!(
            fdbg.node(fdbg.find_node_from_kmer(&kmer(b"nnnT")).unwrap())
                .copy_density(),
            1.0,
        );
    }
    #[test]
    fn em_float_intersection_small() {
        let dbg = mock_intersection_small();
        let genome_size = dbg.genome_size() as CopyDensity;
        let reads = Reads::from(vec![b"ATAGCT".to_vec()]);
        let fdbg = FloatDbg::from_dbg(&dbg);
        let params = PHMMParams::zero_error();
        let (edge_freqs, init_freqs, p) = e_step(&fdbg, &reads, &params);
        m_step(
            &fdbg,
            &edge_freqs,
            &init_freqs,
            genome_size,
            0.1,
            10,
            &params,
        );
    }
}
