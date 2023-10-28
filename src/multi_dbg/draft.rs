//!
//! Constructor of draft dbg from reads or genomes
//!
use super::{CopyNums, MultiDbg};
use crate::common::collection::ReadCollection;
use crate::common::{CopyNum, Freq, Seq};
use crate::distribution::kmer_coverage;
use crate::e2e::Dataset;
use crate::genome::Genome;
use crate::graph::utils::split_node;
use crate::hashdbg::HashDbg;
use crate::hmmv2::{freq::NodeFreqs, hint::Mappings};
use crate::kmer::VecKmer;
use crate::prob::Prob;

use fnv::FnvHashMap as HashMap;
use itertools::{izip, Itertools};
use petgraph::graph::{DiGraph, NodeIndex};
use rustflow::min_flow::{convex::ConvexCost, min_cost_flow_convex_fast, FlowEdge};

///
/// Edge attribute for min_squared_error_copy_nums_from_freqs
///
/// FlowEdge
/// If this node has fixed_copy_num,
/// * demand = fixed_copy_num
/// * capacity = fixed_copy_num
///
/// Otherwise,
/// * demand = 0
/// * capacity = +inf
///
/// ConvexCost
/// * cost = |c - f|^2
///
#[derive(Clone, Debug)]
pub struct MinSquaredErrorCopyNumAndFreq<T: ErrorMetric> {
    ///
    ///
    freqs: Vec<Freq>,
    ///
    ///
    fixed_copy_num: Option<CopyNum>,
    ///
    ///
    non_zero: bool,
    ///
    ///
    phantom: std::marker::PhantomData<T>,
}

impl<T: ErrorMetric> MinSquaredErrorCopyNumAndFreq<T> {
    ///
    /// constructor from Vec<(is_target: bool, freq: Freeq)> and predetermined copy_num
    ///
    pub fn new(freqs: Vec<Freq>, fixed_copy_num: Option<CopyNum>, non_zero: bool) -> Self {
        assert!(
            freqs.iter().all(|f| f.is_finite()),
            "some of freqs in MSE are either nan or infinite"
        );
        MinSquaredErrorCopyNumAndFreq {
            freqs,
            fixed_copy_num,
            non_zero,
            phantom: std::marker::PhantomData,
        }
    }
}

pub trait ErrorMetric: Clone {
    fn cost(freqs: &[Freq], copy_num: CopyNum) -> f64;
}

/// ErrorMetric V1
///
/// h(c) = |c-f|^2
///
#[derive(Clone, Debug)]
pub struct V1Error {}

impl ErrorMetric for V1Error {
    fn cost(freqs: &[Freq], copy_num: CopyNum) -> f64 {
        freqs
            .iter()
            .map(|&freq| (copy_num as f64 - freq).powi(2))
            .sum()
    }
}

/// ErrorMetric V2
///
/// h(c) = |1 - c/f|^2
///
#[derive(Clone, Debug)]
pub struct V2Error {}

impl ErrorMetric for V2Error {
    fn cost(freqs: &[Freq], copy_num: CopyNum) -> f64 {
        freqs
            .iter()
            .map(|&freq| (1.0 - (copy_num as f64 / (freq + 1e-7))).powi(2))
            .sum()
    }
}

/// ErrorMetric V4
///
/// h(c) = |1 - c/f|^2 + |f/c - 1|^2
///
#[derive(Clone, Debug)]
pub struct V4Error {}

impl ErrorMetric for V4Error {
    fn cost(freqs: &[Freq], copy_num: CopyNum) -> f64 {
        freqs
            .iter()
            .map(|&freq| {
                (1.0 - (copy_num as f64 / (freq + 1e-1))).powi(2)
                    + ((freq / (copy_num as f64 + 1e-1)) - 1.0).powi(2)
            })
            .sum()
    }
}

///
/// maximum copy number
///
/// this corresponds to the capacity of edbg min-flow calculation.
///
pub const MAX_COPY_NUM_OF_EDGE: usize = 1000;

impl<T: ErrorMetric> FlowEdge<usize> for MinSquaredErrorCopyNumAndFreq<T> {
    fn demand(&self) -> usize {
        match self.fixed_copy_num {
            Some(fixed_copy_num) => fixed_copy_num,
            None => {
                if self.non_zero {
                    1
                } else {
                    0
                }
            }
        }
    }
    fn capacity(&self) -> usize {
        match self.fixed_copy_num {
            Some(fixed_copy_num) => fixed_copy_num,
            None => MAX_COPY_NUM_OF_EDGE,
        }
    }
}

///
/// Use edbg edge (with a freq) in min-flow.
///
/// if the kmer corresponding to the edge is not emittable, the cost
/// should be ignored.
///
/// `x = copy_num, y = freq`
///
/// # candidate weights
///
/// * `h(x) = |x - f|^2`
/// * `h(x) = |1 - x/(f+e)|^2`
/// * `h(x) = |1 - f/(x+1)|^2` not good
///
impl<T: ErrorMetric> ConvexCost<usize> for MinSquaredErrorCopyNumAndFreq<T> {
    fn convex_cost(&self, copy_num: usize) -> f64 {
        T::cost(&self.freqs, copy_num)
    }
}

///
///
///
#[derive(Clone, Debug, Copy)]
pub enum TerminalCount {
    ///
    ///
    ///
    Free,
    ///
    /// Specify n_haplotypes by splitting terminal node into two nodes
    /// connected by an edge whose demand and capacity is set to the `n_haplotypes`.
    ///
    Fixed(usize),
    ///
    /// Disconnect in-edges and out-edges of terminal node, so that the number of haplotypes will
    /// not be changed by updating with cycle in the residue graph.
    ///
    /// Note that disconnecting terminals breaks flow consistency, so it cannot be used to run
    /// min-flow algorithms. (only for finding cycles and neighbors.)
    ///
    Disconnect,
}

impl MultiDbg {
    ///
    ///
    ///
    pub fn mappings_to_freqs(&self, mappings: &Mappings) -> NodeFreqs {
        // freq is defined for each nodes in hmm (= edges in dbg)
        let mut freqs: NodeFreqs = NodeFreqs::new(self.n_edges_full(), 0.0, true);
        for mapping in mappings {
            for i in 0..mapping.len() {
                for (&node, prob) in izip!(&mapping.nodes[i], &mapping.probs[i]) {
                    freqs[node] += prob.to_value();
                }
            }
        }
        freqs
    }
    ///
    /// * terminal_count:
    /// * not_make_new_zero_edge: if true, demand of non-zero edge will be 1.
    ///
    pub fn to_min_squared_error_copy_nums_network<T: ErrorMetric>(
        &self,
        freqs: &NodeFreqs,
        coverage: f64,
        terminal_count: TerminalCount,
        not_make_new_zero_edge: bool,
    ) -> DiGraph<(), MinSquaredErrorCopyNumAndFreq<T>> {
        let mut net = self.graph_compact().map(
            |_, _| (),
            |edge_in_comapct, _| {
                let freqs_of_edge = self
                    .edges_in_full(edge_in_comapct)
                    .iter()
                    .filter_map(|&edge_in_full| {
                        if self.graph_full()[edge_in_full].is_null_base() {
                            None
                        } else {
                            Some(freqs[NodeIndex::new(edge_in_full.index())] / coverage)
                        }
                    })
                    .collect();
                let copy_num = self.copy_num_of_edge_in_compact(edge_in_comapct);
                let non_zero = not_make_new_zero_edge && copy_num != 0;
                MinSquaredErrorCopyNumAndFreq::new(freqs_of_edge, None, non_zero)
            },
        );

        match terminal_count {
            TerminalCount::Fixed(n_haplotypes) => {
                // split terminal node into two
                let terminal = self
                    .terminal_node_compact()
                    .expect("n_haplotype is specified, but there is no terminal node");
                split_node(
                    &mut net,
                    terminal,
                    Some(MinSquaredErrorCopyNumAndFreq::new(
                        vec![],
                        Some(n_haplotypes),
                        false,
                    )),
                );
            }
            TerminalCount::Disconnect => {
                // split terminal node into two nodes to prohibit cycles passing the terminal node.
                if let Some(terminal) = self.terminal_node_compact() {
                    // split terminal node into two disconnected nodes
                    split_node(&mut net, terminal, None);
                }
            }
            TerminalCount::Free => {
                // no splitting is required.
            }
        }

        // println!("[mse] network");
        // println!("[mse] {:?}", petgraph::dot::Dot::with_config(&net, &[]));
        net
    }
    ///
    /// * n_haplotypes
    ///     if the number of linear haplotypes is known and specified, set to Some(n_haplotypes).
    ///
    pub fn min_squared_error_copy_nums_from_freqs(
        &self,
        freqs: &NodeFreqs,
        coverage: f64,
        n_haplotypes: Option<usize>,
    ) -> CopyNums {
        let terminal_count = if let Some(n_haplotypes) = n_haplotypes {
            TerminalCount::Fixed(n_haplotypes)
        } else {
            TerminalCount::Free
        };
        let net = self.to_min_squared_error_copy_nums_network::<V4Error>(
            freqs,
            coverage,
            terminal_count,
            false,
        );
        let copy_nums = min_cost_flow_convex_fast(&net).expect("mse flownetwork cannot be solved");
        // println!("[mse] copy_nums={}", copy_nums);
        if n_haplotypes.is_some() {
            self.trim_last(&copy_nums)
        } else {
            copy_nums
        }
    }
    ///
    ///
    fn trim_last(&self, copy_nums: &CopyNums) -> CopyNums {
        assert_eq!(copy_nums.len(), self.n_edges_compact() + 1);
        // match size
        let mut ret = CopyNums::new(self.n_edges_compact(), 0);
        for e in self.graph_compact().edge_indices() {
            ret[e] = copy_nums[e];
        }
        ret
    }
}

impl MultiDbg {
    ///
    /// Create from reads V2
    ///
    pub fn create_draft_from_reads_v2<S: Seq>(
        k: usize,
        reads: &ReadCollection<S>,
        p_error: Prob,
        genome_size: usize,
        n_haplotypes: Option<usize>,
        min_count: usize,
        min_deadend_count: usize,
    ) -> Self {
        eprintln!("[draftv2] reads..");
        let mut hd: HashDbg<VecKmer> = HashDbg::from_fragment_seqs(k, reads);
        eprintln!("[draftv2] raw degree_stats={:?}", hd.degree_stats());
        eprintln!("[draftv2] raw count_stats={:?}", hd.copy_num_stats());
        let n_removed = hd.remove_rare_kmers(min_count);
        eprintln!("[draftv2] removed {} rare k-mers", n_removed);
        let n_removed_deadends = hd.remove_deadends(min_deadend_count);
        eprintln!(
            "[draftv2] removed {} rare deadend k-mers",
            n_removed_deadends
        );
        let (starts, ends) = hd.augment_deadends();
        eprintln!(
            "[draftv2] deadends start={} end={}",
            starts.iter().join(","),
            ends.iter().join(","),
        );

        let components = hd.connected_components();
        for (i, component) in components.into_iter().enumerate() {
            eprintln!("component #{} {}", i, component.len());
        }
        let hd = hd.largest_component();

        eprintln!("[draftv2] n={}", hd.n());
        let coverage = reads.coverage(genome_size);
        let adjusted_coverage = kmer_coverage(k, reads.average_length(), coverage, p_error);
        eprintln!(
            "[draftv2] coverage={} adjusted_coverage={}",
            coverage, adjusted_coverage
        );
        let hd =
            hd.generate_hashdbg_with_min_squared_error_copy_nums(adjusted_coverage, n_haplotypes);
        eprintln!("[draftv2] raw degree_stats={:?}", hd.degree_stats());
        eprintln!("[draftv2] raw count_stats={:?}", hd.copy_num_stats());

        Self::from_hashdbg(&hd, false)
    }
    /// Create read-draft MultiDbg from dataset V2
    /// with `min_count=2` and `min_deadend_count=coverage/4`
    pub fn create_draft_from_dataset(k: usize, dataset: &Dataset) -> Self {
        Self::create_draft_from_dataset_with(k, dataset, 2, (dataset.coverage() / 4.0) as usize)
    }
    /// Create read-draft MultiDbg from dataset V2
    /// with specifying `min_count` and `min_deadend_count`
    pub fn create_draft_from_dataset_with(
        k: usize,
        dataset: &Dataset,
        min_count: usize,
        min_deadend_count: usize,
    ) -> Self {
        Self::create_draft_from_reads_v2(
            k,
            dataset.reads(),
            dataset.params().p_error(),
            dataset.genome_size(),
            Some(dataset.genome().n_linear_haplotypes()),
            min_count,
            min_deadend_count,
        )
    }
}

/// assertion
#[allow(dead_code)]
fn check_dbg_contains_true_kmers(dbg: &MultiDbg, genome: &Genome) {
    assert!(dbg.paths_from_styled_seqs(genome).is_ok());
}

#[allow(dead_code)]
fn show_kmer_copy_num_map(map: &HashMap<VecKmer, CopyNum>) {
    for (kmer, copy_num) in map {
        println!("{} {}x", kmer, copy_num);
    }
}

#[allow(dead_code)]
fn show_kmer_copy_num_map_diff(a: &HashMap<VecKmer, CopyNum>, b: &HashMap<VecKmer, CopyNum>) {
    println!("--- diff ---");
    for (x, xa) in a {
        if let Some(xb) = b.get(x) {
            if xa != xb {
                println!("{} a:{}x b:{}x", x, xa, xb);
            }
        } else {
            println!("{} a:{}x b:-", x, xa);
        }
    }
    for (x, xb) in b {
        if !a.contains_key(x) {
            println!("{} a:-  b:{}x", x, xb);
        }
    }
    println!("--- diff ---");
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::e2e::{generate_dataset, ReadType};
    use crate::genome;
    use crate::hmmv2::params::PHMMParams;
    use crate::utils::timer;

    fn check_no_diff_between_draft_v1_and_v2(dataset: Dataset, k: usize) {
        // starts_and_ends_of_genome(&genome, k);
        // let (d1, t1) = timer(|| MultiDbg::create_draft_from_dataset_old(k, &dataset));
        // check_dbg_contains_true_kmers(&d1, dataset.genome());

        let hd: HashDbg<VecKmer> = HashDbg::from_fragment_seqs(k, dataset.reads());

        let (d2, t2) = timer(|| MultiDbg::create_draft_from_dataset(k, &dataset));
        // if missing, dump raw hashdbg
        if !d2.paths_from_styled_seqs(dataset.genome()).is_ok() {
            hd.to_gfa_file(format!("tmp.gfa"));
        }
        check_dbg_contains_true_kmers(&d2, dataset.genome());
        assert_eq!(d2.n_haplotypes(), 2);

        // let m1 = d1.to_kmer_copy_num_map();
        // let m2 = d2.to_kmer_copy_num_map();
        // show_kmer_copy_num_map_diff(&m1, &m2);
        // assert_eq!(d1.to_kmer_set(), d2.to_kmer_set());

        // m1, m2);
        // println!("################## t1={} t2={} ################", t1, t2);
    }

    #[test]
    fn mse_cost() {
        let w = MinSquaredErrorCopyNumAndFreq::<V1Error>::new(vec![], None, false);
        assert_eq!(w.demand(), 0);
        assert_eq!(w.capacity(), MAX_COPY_NUM_OF_EDGE);
        assert_eq!(w.convex_cost(0), 0.0);
        assert_eq!(w.convex_cost(1), 0.0);

        let w = MinSquaredErrorCopyNumAndFreq::<V1Error>::new(vec![1.0], None, false);
        assert_eq!(w.demand(), 0);
        assert_eq!(w.capacity(), MAX_COPY_NUM_OF_EDGE);
        assert_eq!(w.convex_cost(0), 1.0);
        assert_eq!(w.convex_cost(1), 0.0);

        let w = MinSquaredErrorCopyNumAndFreq::<V1Error>::new(vec![1.0, 2.0], None, false);
        assert_eq!(w.demand(), 0);
        assert_eq!(w.capacity(), MAX_COPY_NUM_OF_EDGE);
        assert_eq!(w.convex_cost(0), 1.0 + 4.0);
        assert_eq!(w.convex_cost(1), 0.0 + 1.0);

        let w = MinSquaredErrorCopyNumAndFreq::<V1Error>::new(vec![1.0], Some(2), false);
        assert_eq!(w.demand(), 2);
        assert_eq!(w.capacity(), 2);
        assert_eq!(w.convex_cost(2), 1.0);
    }
    #[test]
    #[ignore]
    fn draft_check_v2_n_1_to_12() {
        for n in 1..=12 {
            println!("################## n={} ################", n);
            let genome = genome::u500(n);
            let k = 40;
            let dataset = generate_dataset(
                genome.clone(),
                0,
                20,
                1000,
                ReadType::FragmentWithRevComp,
                PHMMParams::uniform(0.001),
            );
            check_no_diff_between_draft_v1_and_v2(dataset, k);
        }
    }
    #[test]
    fn draft_check_v2_n6() {
        let genome = genome::u500(6);
        let k = 40;
        let dataset = generate_dataset(
            genome.clone(),
            0,
            20,
            1000,
            ReadType::FragmentWithRevComp,
            PHMMParams::uniform(0.001),
        );
        check_no_diff_between_draft_v1_and_v2(dataset, k);
    }
}
