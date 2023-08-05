//!
//! Constructor of draft dbg from reads or genomes
//!
use super::{CopyNums, MultiDbg};
use crate::common::collection::starts_and_ends_of_genome;
use crate::common::{CopyNum, Freq, Seq, StyledSequence};
use crate::dbg::{draft::EndNodeInference, Dbg, SimpleDbg};
use crate::e2e::Dataset;
use crate::graph::utils::split_node;
use crate::hmmv2::{freq::NodeFreqs, hint::Mappings};
use crate::kmer::VecKmer;
use crate::utils::timer;
use itertools::izip;
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
pub struct MinSquaredErrorCopyNumAndFreq {
    ///
    ///
    freqs: Vec<Freq>,
    ///
    ///
    fixed_copy_num: Option<CopyNum>,
    ///
    ///
    non_zero: bool,
}

impl MinSquaredErrorCopyNumAndFreq {
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
        }
    }
}

///
/// maximum copy number
///
/// this corresponds to the capacity of edbg min-flow calculation.
///
pub const MAX_COPY_NUM_OF_EDGE: usize = 1000;

impl FlowEdge<usize> for MinSquaredErrorCopyNumAndFreq {
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
impl ConvexCost<usize> for MinSquaredErrorCopyNumAndFreq {
    fn convex_cost(&self, copy_num: usize) -> f64 {
        self.freqs
            .iter()
            // .map(|&freq| (copy_num as f64 - freq).powi(2))
            .map(|&freq| (1.0 - (copy_num as f64 / (freq + 0.0000001))).powi(2))
            .sum()
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
    pub fn to_min_squared_error_copy_nums_network(
        &self,
        freqs: &NodeFreqs,
        coverage: f64,
        terminal_count: TerminalCount,
        not_make_new_zero_edge: bool,
    ) -> DiGraph<(), MinSquaredErrorCopyNumAndFreq> {
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
        let net =
            self.to_min_squared_error_copy_nums_network(freqs, coverage, terminal_count, false);
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
    /// Create from reads
    ///
    pub fn create_draft_from_reads<T>(
        k: usize,
        seqs: T,
        base_coverage: f64,
        ave_read_length: usize,
        p_error: f64,
        end_node_inference: &EndNodeInference<VecKmer>,
    ) -> Self
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        let dbg: SimpleDbg<VecKmer> =
            SimpleDbg::create_draft_from_fragment_seqs_with_adjusted_coverage(
                k,
                seqs,
                base_coverage,
                ave_read_length,
                p_error,
                end_node_inference,
            );
        dbg.into()
    }
    ///
    /// Create from styled seqs
    ///
    pub fn create_from_styled_seqs<T>(k: usize, seqs: T) -> Self
    where
        T: IntoIterator,
        T::Item: AsRef<StyledSequence>,
    {
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_styled_seqs(k, seqs);
        dbg.into()
    }
    ///
    /// Create read-draft MultiDbg from dataset
    ///
    /// * use end inference from genome
    /// * true path
    ///
    pub fn create_draft_from_dataset(k: usize, dataset: &Dataset) -> Self {
        let d = Self::create_draft_from_reads(
            k,
            dataset.reads(),
            dataset.coverage(),
            dataset.average_read_length(),
            dataset.params().p_error().to_value(),
            // &EndNodeInference::Auto,
            &EndNodeInference::Custom(starts_and_ends_of_genome(dataset.genome(), k)),
        );
        d
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome;

    #[test]
    #[ignore]
    fn from_styled_seqs_large() {
        let genome = genome::tandem_repeat_polyploid_with_unique_homo_ends(
            1_000, 1_000, 0, 0.0, 0, 1_000, 2, 0.01, 0,
        );
        let (mdbg, t) = timer(|| MultiDbg::create_from_styled_seqs(40, &genome));
        // ~2391ms
        println!("created mdbg in {}ms", t);
        mdbg.to_gfa_file("g1m.gfa");
        let (m, t) = timer(|| mdbg.to_kmer_map());
        // ~182ms
        println!("created map in {}ms", t);
        println!("m {}", m.len());

        let (n, t) = timer(|| mdbg.n_euler_circuits());
        println!("n_euler={} in {}ms", n, t);
    }
    #[test]
    fn mse_cost() {
        let w = MinSquaredErrorCopyNumAndFreq::new(vec![], None, false);
        assert_eq!(w.demand(), 0);
        assert_eq!(w.capacity(), MAX_COPY_NUM_OF_EDGE);
        assert_eq!(w.convex_cost(0), 0.0);
        assert_eq!(w.convex_cost(1), 0.0);

        let w = MinSquaredErrorCopyNumAndFreq::new(vec![1.0], None, false);
        assert_eq!(w.demand(), 0);
        assert_eq!(w.capacity(), MAX_COPY_NUM_OF_EDGE);
        assert_eq!(w.convex_cost(0), 1.0);
        assert_eq!(w.convex_cost(1), 0.0);

        let w = MinSquaredErrorCopyNumAndFreq::new(vec![1.0, 2.0], None, false);
        assert_eq!(w.demand(), 0);
        assert_eq!(w.capacity(), MAX_COPY_NUM_OF_EDGE);
        assert_eq!(w.convex_cost(0), 1.0 + 4.0);
        assert_eq!(w.convex_cost(1), 0.0 + 1.0);

        let w = MinSquaredErrorCopyNumAndFreq::new(vec![1.0], Some(2), false);
        assert_eq!(w.demand(), 2);
        assert_eq!(w.capacity(), 2);
        assert_eq!(w.convex_cost(2), 1.0);
    }
}
