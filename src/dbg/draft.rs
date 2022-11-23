//!
//! Constructor of draft dbg
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::common::{ei, ni, CopyNum, Freq, Seq};
use crate::dbg::edge_centric::compact::compacted_flow_into_original_flow;
use crate::hmmv2::freq::NodeFreqs;
use crate::kmer::kmer::KmerLike;
use crate::min_flow::{
    convex::ConvexCost, flow::FlowEdge, min_cost_flow_convex_fast, total_cost, Cost,
};

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    /// Create draft dbg from reads
    ///
    /// 1. remove 1x-copy nodes
    /// 2. remove deadend nodes
    /// 3. assign approximated (flow-constraint-satisfied) copy_nums
    ///
    pub fn create_draft_from_seqs<T>(k: usize, seqs: T, coverage: f64) -> Self
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        eprintln!("[draft] constructing base");
        let mut dbg = Self::from_seqs(k, seqs);
        eprintln!("[draft] n_nodes0={}", dbg.n_nodes());
        eprintln!("[draft] n_edges0={}", dbg.n_edges());
        // 1&2
        dbg.remove_nodes(2);
        dbg.remove_deadend_nodes();
        eprintln!("[draft] n_nodes={}", dbg.n_nodes());
        eprintln!("[draft] n_edges={}", dbg.n_edges());
        eprintln!("[draft] copy_num_stats_raw={:?}", dbg.copy_num_stats());
        eprintln!("[draft] degree_stats={:?}", dbg.degree_stats());
        // 3
        let freqs = dbg.to_node_freqs() / coverage as f64;
        dbg.set_copy_nums_all_zero();
        let (copy_nums_approx, cost) = dbg
            .min_squared_error_copy_nums_from_freqs_compacted(&freqs)
            .unwrap();
        eprintln!("[draft] approx_cost={}", cost);
        dbg.set_node_copy_nums(&copy_nums_approx);
        eprintln!("[draft] copy_num_stats_approx={:?}", dbg.copy_num_stats());
        dbg
    }
    /// Create draft dbg from fragment reads
    ///
    /// WIP
    ///
    pub fn create_draft_from_fragment_seqs<T>(k: usize, seqs: T, coverage: f64) -> Self
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        eprintln!("[draft_frag] constructing raw dbg");
        let mut dbg = Self::from_fragment_seqs(k, seqs);
        assert!(dbg.has_no_duplicated_node());
        eprintln!("[draft_frag] n_nodes_raw={}", dbg.n_nodes());
        eprintln!("[draft_frag] n_edges_raw={}", dbg.n_edges());
        eprintln!("[draft_frag] copy_num_stats_raw={:?}", dbg.copy_num_stats());
        eprintln!("[draft_frag] degree_stats_raw={:?}", dbg.degree_stats());
        dbg.remove_nodes(2);
        eprintln!("[draft_frag] n_nodes_raw={}", dbg.n_nodes());
        eprintln!("[draft_frag] n_edges_raw={}", dbg.n_edges());
        dbg
    }
    ///
    /// by solving min-cost-flow
    ///
    pub fn min_squared_error_copy_nums_from_freqs(
        &self,
        freqs: &NodeFreqs,
    ) -> Option<(NodeCopyNums, Cost)> {
        let graph = self.to_edbg_graph(
            |_| (),
            |v, weight| MinSquaredErrorCopyNumAndFreq::new(vec![(weight.is_emittable(), freqs[v])]),
        );
        min_cost_flow_convex_fast(&graph).map(|flow| {
            let cost = total_cost(&graph, &flow);
            // an edge in edbg corresponds to a node in dbg
            // so edgevec for edbg can be converted to nodevec for dbg.
            (flow.switch_index(), cost)
        })
    }
    ///
    /// by solving min-cost-flow
    /// uses compacted edbg
    ///
    pub fn min_squared_error_copy_nums_from_freqs_compacted(
        &self,
        freqs: &NodeFreqs,
    ) -> Option<(NodeCopyNums, Cost)> {
        let graph = self.to_compact_edbg_graph();
        // println!("{}", petgraph::dot::Dot::with_config(&graph, &[]));
        let flow_network = graph.map(
            |_, _| (),
            |_, weight| {
                let freqs = weight
                    .origin_edges()
                    .iter()
                    .map(|&edge| {
                        let node = ni(edge.index());
                        (self.is_emittable(node), freqs[node])
                    })
                    .collect();
                MinSquaredErrorCopyNumAndFreq::new(freqs)
            },
        );
        min_cost_flow_convex_fast(&flow_network).map(|flow| {
            let cost = total_cost(&flow_network, &flow);
            // an edge in edbg corresponds to a node in dbg
            // so edgevec for edbg can be converted to nodevec for dbg.
            // convert flow in compacted-edbg into flow in original edbg
            let flow_in_original = compacted_flow_into_original_flow(self.n_nodes(), &graph, &flow);
            (flow_in_original.switch_index(), cost)
        })
    }
    ///
    /// get current copy nums as node freqs
    ///
    pub fn to_node_freqs(&self) -> NodeFreqs {
        let copy_nums = self.to_node_copy_nums();
        NodeFreqs::from_inner_vec(
            copy_nums
                .to_inner_vec()
                .into_iter()
                .map(|copy_num| copy_num as Freq)
                .collect(),
        )
    }
}

///
/// Edge attribute for min_squared_error_copy_nums_from_freqs
///
/// FlowEdge
/// * demand = 0
/// * capacity = +inf
///
/// ConvexCost
/// * cost = |c - f|^2
///
#[derive(Clone, Debug)]
struct MinSquaredErrorCopyNumAndFreq {
    freqs: Vec<(bool, Freq)>,
}

impl MinSquaredErrorCopyNumAndFreq {
    ///
    /// constructor
    ///
    pub fn new(freqs: Vec<(bool, Freq)>) -> Self {
        MinSquaredErrorCopyNumAndFreq { freqs }
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
        0
    }
    fn capacity(&self) -> usize {
        MAX_COPY_NUM_OF_EDGE
    }
}

///
/// Use edbg edge (with a freq) in min-flow.
///
/// if the kmer corresponding to the edge is not emittable, the cost
/// should be ignored.
///
impl ConvexCost<usize> for MinSquaredErrorCopyNumAndFreq {
    fn convex_cost(&self, copy_num: usize) -> f64 {
        self.freqs
            .iter()
            .map(|&(is_target, freq)| {
                if is_target {
                    (copy_num as f64 - freq).powi(2)
                } else {
                    0.0
                }
            })
            .sum()
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks;
    use crate::dbg::SimpleDbg;
    use crate::e2e::{generate_simple_genome_fragment_dataset, generate_simple_genome_mock};
    use crate::io::write_string;
    use crate::kmer::VecKmer;

    #[test]
    fn dbg_to_freqs_test() {
        let dbg = mocks::mock_intersection_small();
        let nc = dbg.to_node_copy_nums();
        let nf = dbg.to_node_freqs();
        println!("nc={}", nc);
        println!("nf={}", nf);
        assert_eq!(nc.to_vec(), vec![1; dbg.n_nodes()]);
        assert_eq!(nf.to_vec(), vec![1.0; dbg.n_nodes()]);
    }
    #[test]
    fn dbg_min_squared_error_copy_nums_obvious_test() {
        let dbg = mocks::mock_intersection();
        let nc = dbg.to_node_copy_nums();
        let nf = dbg.to_node_freqs();
        // case 1: true freq
        let (fitted, cost) = dbg.min_squared_error_copy_nums_from_freqs(&nf).unwrap();
        assert_eq!(fitted, nc);
        assert_eq!(cost, 0.0);
        // case 2: all zero freq
        let nf0 = NodeFreqs::new(dbg.n_nodes(), 0.0);
        let (fitted, cost) = dbg.min_squared_error_copy_nums_from_freqs(&nf0).unwrap();
        assert_eq!(fitted.sum(), 0);
        assert_eq!(cost, 0.0);
    }
    #[test]
    fn dbg_min_squared_error_copy_nums_simple_genome_test() {
        let experiment = generate_simple_genome_mock();
        let dbg_raw = experiment.dbg_raw.clone();
        // freq = (read occurrences of kmers) / (coverage)
        let freq = dbg_raw.to_node_freqs() / 20.0;
        // TODO show histogram of coverages
        let (copy_nums_true, _) = dbg_raw
            .to_copy_nums_of_styled_seqs(experiment.genome())
            .unwrap();

        // (1) normal
        let (approx, cost) = dbg_raw
            .min_squared_error_copy_nums_from_freqs(&freq)
            .unwrap();
        println!("{}", freq);
        println!("{}", approx);
        println!("{}", copy_nums_true);
        println!("{}", approx.dist(&copy_nums_true));
        println!("{}", cost);
        assert_eq!(approx.dist(&copy_nums_true), 0);
        assert_eq!(approx, copy_nums_true);
        assert!(4.7 <= cost && cost <= 4.8);

        // (2) compact
        let (approx2, cost2) = dbg_raw
            .min_squared_error_copy_nums_from_freqs_compacted(&freq)
            .unwrap();
        assert_eq!(approx2.dist(&copy_nums_true), 0);
        assert_eq!(approx2, approx);
        assert!(4.7 <= cost2 && cost2 <= 4.8);
    }
    #[test]
    fn dbg_create_draft_simple_genome_test() {
        let experiment = generate_simple_genome_mock();
        let copy_nums_true = experiment.dbg_draft_true.unwrap().to_node_copy_nums();
        let approx = experiment.dbg_draft.unwrap().to_node_copy_nums();
        assert_eq!(approx.dist(&copy_nums_true), 0);
    }
    #[test]
    fn dbg_create_draft_fragment_test() {
        let dataset = generate_simple_genome_fragment_dataset();
        dataset.show_genome();
        dataset.show_reads();
        println!("coverage={}", dataset.coverage());
        let dbg: SimpleDbg<VecKmer> =
            SimpleDbg::create_draft_from_fragment_seqs(32, dataset.reads(), dataset.coverage());
        let json = dbg.to_cytoscape();
        write_string("draft_from_fragment.json", &json).unwrap();
    }
}
