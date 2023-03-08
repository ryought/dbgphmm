//!
//! Compression V3
//! By Naive main factor improving
//!
use super::kmer_info::{create_kmer_infos as create_plain_kmer_infos, KmerInfo};
use super::q::{q_score_clamped, QScore};
use crate::common::{CopyNum, Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::e2e::Experiment;
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::{EdgeFreqs, NodeFreqs};
use crate::kmer::kmer::KmerLike;
use crate::min_flow::residue::{
    improve_flow_convex_with_update_info, CycleDetectMethod, ResidueDirection, UpdateInfo,
};
use crate::min_flow::utils::clamped_log_with;
use crate::min_flow::{min_cost_flow_from_convex, total_cost, Cost};
use crate::prob::Prob;
use crate::vector::{DenseStorage, EdgeVec, NodeVec, Storage};
use petgraph::graph::{EdgeIndex, NodeIndex};

use crate::dbg::draft::MAX_COPY_NUM_OF_EDGE;
use crate::dbg::edge_centric::impls::SimpleEDbgEdgeWithAttr;
use crate::min_flow::convex::ConvexCost;
use crate::min_flow::FlowEdge;

// use same e-step function in V2
pub use crate::em::compression::v2::e_step;

#[derive(Clone, Debug, Copy, PartialEq, Default)]
pub struct CompressionV3KmerInfo {
    kmer_info: KmerInfo,
    zero_penalty: f64,
}

pub type V3KmerInfos = NodeVec<DenseStorage<CompressionV3KmerInfo>>;
pub type SimpleEDbgEdgeWithV3KmerInfos<K> = SimpleEDbgEdgeWithAttr<K, CompressionV3KmerInfo>;

#[derive(Clone, Debug, Copy, PartialEq)]
pub enum Diff {
    /// +1
    Inc,
    /// -1
    Dec,
    /// 0
    Nop,
}

impl Diff {
    ///
    /// Apply the diff operation
    ///
    pub fn apply(&self, value: usize) -> usize {
        match self {
            Diff::Inc => value.saturating_add(1),
            Diff::Dec => value.saturating_sub(1),
            Diff::Nop => value,
        }
    }
}

impl CompressionV3KmerInfo {
    ///
    /// Current copy number of the kmer
    ///
    pub fn copy_num(&self) -> CopyNum {
        self.kmer_info.copy_num
    }
    ///
    /// Exact Q function score
    ///
    pub fn score_from_copy_num(&self, copy_num: CopyNum) -> f64 {
        let diff = if copy_num == self.copy_num() {
            Diff::Nop
        } else if copy_num > self.copy_num() {
            Diff::Inc
        } else {
            Diff::Dec
        };
        self.score(diff)
    }
    ///
    /// Exact Q function score
    ///
    pub fn score(&self, diff: Diff) -> f64 {
        self.x(diff) + self.y(diff) + self.z(diff) + self.r(diff)
    }
    ///
    /// `X_l log (copy_num)` term
    ///
    pub fn x(&self, diff: Diff) -> f64 {
        let new_copy_num = diff.apply(self.kmer_info.copy_num);
        -self.kmer_info.freq * clamped_log_with(new_copy_num, self.zero_penalty)
    }
    ///
    /// `Y_i log (copy_num_intersection)` term
    ///
    pub fn y(&self, diff: Diff) -> f64 {
        let new_copy_num_intersection = diff.apply(self.kmer_info.copy_num_intersection);
        self.kmer_info.freq_intersection
            * clamped_log_with(new_copy_num_intersection, self.zero_penalty)
    }
    ///
    /// `Z log (copy_num_total)` term
    ///
    pub fn z(&self, diff: Diff) -> f64 {
        let new_copy_num_total = diff.apply(self.kmer_info.copy_num_total);
        self.kmer_info.freq_init * clamped_log_with(new_copy_num_total, self.zero_penalty)
    }
    ///
    /// regularization term `R`
    ///
    pub fn r(&self, diff: Diff) -> f64 {
        let new_copy_num_total = diff.apply(self.kmer_info.copy_num_total);
        self.kmer_info.penalty_weight
            * new_copy_num_total
                .abs_diff(self.kmer_info.copy_num_total_expected)
                .pow(2) as f64
    }
}

impl CompressionV3KmerInfo {
    fn new(kmer_info: KmerInfo, zero_penalty: f64) -> Self {
        CompressionV3KmerInfo {
            kmer_info,
            zero_penalty,
        }
    }
}

impl std::fmt::Display for CompressionV3KmerInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "+1={:.4}({:.4},{:.4},{:.4},{:.4}) -1={:.4}({:.4},{:.4},{:.4},{:.4}) ({})",
            // +1
            self.score(Diff::Inc) - self.score(Diff::Nop),
            self.x(Diff::Inc) - self.x(Diff::Nop),
            self.y(Diff::Inc) - self.y(Diff::Nop),
            self.z(Diff::Inc) - self.z(Diff::Nop),
            self.r(Diff::Inc) - self.r(Diff::Nop),
            // -1
            self.score(Diff::Dec) - self.score(Diff::Nop),
            self.x(Diff::Dec) - self.x(Diff::Nop),
            self.y(Diff::Dec) - self.y(Diff::Nop),
            self.z(Diff::Dec) - self.z(Diff::Nop),
            self.r(Diff::Dec) - self.r(Diff::Nop),
            self.kmer_info,
        )
    }
}

impl<K: KmerLike> FlowEdge<usize> for SimpleEDbgEdgeWithV3KmerInfos<K> {
    ///
    /// demand is set to be current copy_num - 1
    ///
    fn demand(&self) -> usize {
        if self.attribute.copy_num() > 0 {
            self.attribute.copy_num() - 1
        } else {
            0
        }
    }
    ///
    /// capacity is set to be current copy_num + 1
    ///
    fn capacity(&self) -> usize {
        if self.attribute.copy_num() < MAX_COPY_NUM_OF_EDGE {
            self.attribute.copy_num() + 1
        } else {
            MAX_COPY_NUM_OF_EDGE
        }
    }
}

impl<K: KmerLike> ConvexCost<usize> for SimpleEDbgEdgeWithV3KmerInfos<K> {
    ///
    ///
    ///
    fn convex_cost(&self, flow: usize) -> f64 {
        self.attribute.score_from_copy_num(flow)
    }
}

///
/// create plain KmerInfos and convert it to V3KmerInfos
///
fn create_kmer_infos<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
    zero_penalty: f64,
) -> V3KmerInfos {
    let mut ki = V3KmerInfos::new(dbg.n_nodes(), CompressionV3KmerInfo::default());
    let ki0 = create_plain_kmer_infos(dbg, edge_freqs, init_freqs, genome_size, penalty_weight);
    for (node, _) in dbg.nodes() {
        ki[node] = CompressionV3KmerInfo::new(ki0[node], zero_penalty);
    }
    ki
}

///
/// result of m-step oneshot run
///
/// * Update
/// * NoImprove
/// * NoNegCycle
///
pub enum MStepResult<N: DbgNode, E: DbgEdge> {
    /// The found negative cycle improved q-score
    Update(Dbg<N, E>, QScore, Cost, Updates),
    /// A negative cycle was found, but it did not improve q-score
    NoImprove(QScore, Cost, Updates),
    /// Any negative cycle was not found.
    NoNegCycle,
}

impl<N: DbgNode, E: DbgEdge> MStepResult<N, E> {
    ///
    /// is MStepResult::Update or not?
    ///
    pub fn is_update(&self) -> bool {
        match self {
            MStepResult::Update(_, _, _, _) => true,
            _ => false,
        }
    }
    ///
    /// return contents of MStepResult::Update
    ///
    pub fn unwrap(self) -> (Dbg<N, E>, QScore, Cost, Updates) {
        match self {
            MStepResult::Update(d, q, c, u) => (d, q, c, u),
            _ => panic!(),
        }
    }
}

pub type Updates = Vec<(NodeIndex, ResidueDirection)>;

///
/// conversion from UpdateInfo (= edge index of edbg and its up/down)
/// into node UpdateInfo
///
fn edbg_update_info_to_kmer_update_info(
    update_info: &[(EdgeIndex, ResidueDirection)],
) -> Vec<(NodeIndex, ResidueDirection)> {
    update_info
        .iter()
        .map(|(e, dir)| (NodeIndex::new(e.index()), *dir))
        .collect()
}

///
/// Run a single update
///
fn m_step_once<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
    zero_penalty: f64,
) -> MStepResult<N, E> {
    // (0) calculate original q_score
    let q_score = q_score_clamped(
        &dbg,
        &edge_freqs,
        &init_freqs,
        genome_size,
        penalty_weight,
        zero_penalty,
    );

    // (1) construct edbg with KmerInfo
    let infos = create_kmer_infos(
        dbg,
        edge_freqs,
        init_freqs,
        genome_size,
        penalty_weight,
        zero_penalty,
    );
    let edbg = dbg.to_edbg_with_attr(Some(&infos));
    // dbg.draw_with_vecs(&[&infos], &[]);
    // dbg.draw_plain_with_vecs(&[&infos], &[]);

    // (2) min-flow optimization starts from current copy nums
    // this edbg is not convex, so use MMWC
    let original_copy_nums = dbg.to_node_copy_nums().switch_index();
    match improve_flow_convex_with_update_info(
        &edbg.graph,
        &original_copy_nums,
        CycleDetectMethod::MinMeanWeightCycle,
    ) {
        Some((copy_nums, update_info)) => {
            // approximated cost difference
            let cost_diff =
                total_cost(&edbg.graph, &copy_nums) - total_cost(&edbg.graph, &original_copy_nums);

            // collect copy number updates
            let updates = edbg_update_info_to_kmer_update_info(&update_info);

            // calculate actual q-score
            let mut dbg_new = dbg.clone();
            dbg_new.set_node_copy_nums(&copy_nums.switch_index());
            let q_score_new = q_score_clamped(
                &dbg_new,
                &edge_freqs,
                &init_freqs,
                genome_size,
                penalty_weight,
                zero_penalty,
            );

            // if new q_score is bigger, accept the change.
            if q_score_new.total() > q_score.total() {
                MStepResult::Update(dbg_new, q_score_new, cost_diff, updates)
            } else {
                MStepResult::NoImprove(q_score_new, cost_diff, updates)
            }
        }
        None => MStepResult::NoNegCycle,
    }
}

fn show_updates<N: DbgNode, E: DbgEdge>(dbg: &Dbg<N, E>, updates: &Updates) {
    for (node, direction) in updates {
        match direction {
            ResidueDirection::Up => {
                println!("+1 {}", dbg.kmer(*node));
            }
            ResidueDirection::Down => {
                println!("-1 {}", dbg.kmer(*node));
            }
        }
    }
}

///
/// M-step of compression_v3
///
/// find a candidate update of current flow
/// using variational approximation
///
fn m_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    edge_freqs: &EdgeFreqs,
    init_freqs: &NodeFreqs,
    genome_size: CopyNum,
    penalty_weight: f64,
    zero_penalty: f64,
    n_max_update: usize,
) -> Vec<(Dbg<N, E>, QScore, Cost)> {
    let mut dbg_current = dbg.clone();
    let mut ret = Vec::new();

    for i in 0..n_max_update {
        // try to improve
        match m_step_once(
            &dbg_current,
            edge_freqs,
            init_freqs,
            genome_size,
            penalty_weight,
            zero_penalty,
        ) {
            MStepResult::Update(dbg_new, q_score, cost, updates) => {
                // found better dbg!
                // show_updates(&dbg_current, &updates);
                ret.push((dbg_new.clone(), q_score, cost));
                dbg_current = dbg_new;
            }
            MStepResult::NoImprove(_, _, _) => {
                // not found
                break;
            }
            MStepResult::NoNegCycle => {
                // not found
                break;
            }
        }
    }

    ret
}

///
/// Do compression v3
///
pub fn compression_step<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    genome_size: CopyNum,
    penalty_weight: f64,
    zero_penalty: f64,
    n_max_update: usize,
) -> (Dbg<N, E>, bool, CompressionV3Log<N, E>) {
    // [1] e-step
    let (edge_freqs, init_freqs, p) = e_step(dbg, reads, params);

    // calculate current q-score
    let q_score = q_score_clamped(
        dbg,
        &edge_freqs,
        &init_freqs,
        genome_size,
        penalty_weight,
        zero_penalty,
    );

    // [2] m-step
    let mut steps = m_step(
        dbg,
        &edge_freqs,
        &init_freqs,
        genome_size,
        penalty_weight,
        zero_penalty,
        n_max_update,
    );

    // [3] history
    if steps.len() > 0 {
        let (mut dbg_updated, q_score_updated, cost) = steps.pop().unwrap();
        dbg_updated.remove_zero_copy_node();
        let log = CompressionV3Log::new(
            p,
            q_score,
            q_score_updated,
            cost,
            true,
            dbg_updated.clone(),
            penalty_weight,
            zero_penalty,
        );
        (dbg_updated, true, log)
    } else {
        // dbg is not updated
        let log = CompressionV3Log::new(
            p,
            q_score,
            q_score,
            0.0, //TODO
            false,
            dbg.clone(),
            penalty_weight,
            zero_penalty,
        );
        (dbg.clone(), false, log)
    }
}

///
/// Compression full algorithm by running `compression_step` iteratively.
///
/// * max_iter: max iteration loop count of EM.
/// * on_iter: callback of each end of iteration
///
pub fn compression<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    genome_size: CopyNum,
    penalty_weight: f64,
    zero_penalty: f64,
    n_max_update: usize,
    max_iter: usize,
) -> (Dbg<N, E>, Vec<CompressionV3Log<N, E>>) {
    compression_with_callback(
        dbg,
        reads,
        params,
        genome_size,
        penalty_weight,
        zero_penalty,
        n_max_update,
        max_iter,
        |_, _, _| {},
    )
}

///
/// Compression full algorithm by running `compression_step` iteratively.
///
/// * max_iter: max iteration loop count of EM.
/// * on_iter: callback of each end of iteration
///
pub fn compression_with_callback<
    N: DbgNode,
    E: DbgEdge,
    F: Fn(usize, &Dbg<N, E>, &CompressionV3Log<N, E>),
>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    genome_size: CopyNum,
    penalty_weight: f64,
    zero_penalty: f64,
    n_max_update: usize,
    max_iter: usize,
    on_iter: F,
) -> (Dbg<N, E>, Vec<CompressionV3Log<N, E>>) {
    let mut dbg = dbg.clone();
    let mut logs = Vec::new();

    // iterate EM steps
    for i in 0..max_iter {
        let (dbg_new, is_updated, log) = compression_step(
            &dbg,
            reads,
            params,
            genome_size,
            penalty_weight,
            zero_penalty,
            n_max_update,
        );
        // run callback
        on_iter(i, &dbg_new, &log);

        // save the log
        logs.push(log);

        // if the single EM step does not change the DBG model, stop iteration.
        if !is_updated {
            break;
        }
        dbg = dbg_new;
    }

    (dbg, logs)
}

///
/// Log information store of each iteration in compression v3
///
#[derive(Clone)]
pub struct CompressionV3Log<N: DbgNode, E: DbgEdge> {
    /// Full probability P(R|G) in e-step
    /// that is P(R|G) before improving q-score-function
    pub p: Prob,
    /// q-score before m-step
    pub q0: QScore,
    /// q-score after m-step
    pub q1: QScore,
    /// cost improvement with variational-approximated q-score
    pub cost_diff: Cost,
    /// the dbg changed by this compression step?
    pub is_updated: bool,
    /// resulting dbg
    pub dbg: Dbg<N, E>,
    /// parameter backup of compressionv3
    /// penalty weight (lambda)
    pub penalty_weight: f64,
    /// parameter backup of compressionv3
    /// zero penalty (L0)
    pub zero_penalty: f64,
}

impl<N: DbgNode, E: DbgEdge> CompressionV3Log<N, E> {
    pub fn new(
        p: Prob,
        q0: QScore,
        q1: QScore,
        cost_diff: Cost,
        is_updated: bool,
        dbg: Dbg<N, E>,
        penalty_weight: f64,
        zero_penalty: f64,
    ) -> Self {
        CompressionV3Log {
            p,
            q0,
            q1,
            cost_diff,
            is_updated,
            dbg,
            penalty_weight,
            zero_penalty,
        }
    }
}

//
// to_string of CompressionV3Log
//

impl<N: DbgNode, E: DbgEdge> std::fmt::Display for CompressionV3Log<N, E> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\tis_updated={}\t{}",
            self.p, self.q0, self.q1, self.cost_diff, self.is_updated, self.dbg
        )
    }
}

impl<N: DbgNode, E: DbgEdge> CompressionV3Log<N, E> {
    ///
    /// Detailed output when the origin datset (especially the true genome)
    /// data is available.
    ///
    pub fn to_benchmark_string(&self, dataset: &Experiment) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.dbg
                .benchmark_compression(&dataset, self.penalty_weight),
            self.p.to_log_value(),
            self.q0,
            self.q1,
            self.cost_diff,
            self.dbg,
        )
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;
    use crate::dbg::mocks::*;
    use crate::min_flow::utils::DEFAULT_CLAMP_VALUE;

    #[test]
    fn em_compression_v3_m() {
        let dbg = mock_intersection();
        let reads = Reads {
            reads: vec![
                b"AACTAGCTT".to_vec(),
                b"AACTAGCTT".to_vec(),
                b"AACTAGCTT".to_vec(),
            ],
        };
        let params = PHMMParams::default();
        println!("{}", dbg);
        println!("{}", dbg.n_traverse_choices());
        let (ef, nf, p) = e_step(&dbg, &reads, &params);
        println!("{}", ef);
        println!("{}", nf);

        let lambda = 0.0;
        let genome_size = 9;

        let infos = create_kmer_infos(&dbg, &ef, &nf, genome_size, lambda, DEFAULT_CLAMP_VALUE);
        // dbg.draw_with_vecs(&[&nf], &[&ef]);
        dbg.draw_with_vecs(&[&infos], &[]);
        // println!("{}", infos);
        // let qs = q_score(&dbg, &ef, &nf, genome_size, lambda);
        // println!("{:?}", qs);

        let r = m_step_once(&dbg, &ef, &nf, genome_size, lambda, DEFAULT_CLAMP_VALUE);
        assert!(r.is_update());
        let (dbg, q_score, cost, _) = r.unwrap();
        let copy_nums = dbg.to_node_copy_nums();
        println!("{}", copy_nums);
        println!("{}", cost);
        assert_eq!(
            copy_nums.to_vec(),
            vec![0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0]
        );
    }

    #[test]
    fn em_compression_v3_step() {
        let dbg = mock_intersection();
        let reads = Reads {
            reads: vec![
                b"AACTAGCTT".to_vec(),
                b"AACTAGCTT".to_vec(),
                b"AACTAGCTT".to_vec(),
            ],
        };
        let params = PHMMParams::default();
        println!("dbg0={}", dbg);

        let (new_dbg, is_updated, log) =
            compression_step(&dbg, &reads, &params, 9, 0.0, DEFAULT_CLAMP_VALUE, 1);
        println!("dbg1={}", new_dbg);
    }
}
