//!
//! Constructor of draft dbg
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::common::{ei, ni, CopyNum, Freq, Seq};
use crate::dbg::edge_centric::compact::compacted_flow_into_original_flow;
use crate::distribution::kmer_coverage;
use crate::hmmv2::freq::NodeFreqs;
use crate::kmer::kmer::KmerLike;
use crate::min_flow::{convex::ConvexCost, min_cost_flow_convex_fast, total_cost, Cost, FlowEdge};
use crate::multi_dbg::draft::{MinSquaredErrorCopyNumAndFreq, V2Error};
use crate::prob::Prob;
use crate::utils::timer;
use crate::vector::graph::flow_to_edgevec;
use fnv::FnvHashMap as HashMap;
use petgraph::graph::NodeIndex;

///
///
///
pub enum EndNodeInference<K: KmerLike> {
    ///
    /// End node (source and sink) is automatically infered from reads
    ///
    Auto,
    ///
    /// specify the end nodes
    ///
    Custom((Vec<(K, CopyNum)>, Vec<(K, CopyNum)>)),
}

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
        let (_, t) = timer(|| {
            dbg.remove_nodes(2);
        });
        eprintln!("[draft/remove_1x_nodes] {}", t);
        let (_, t) = timer(|| {
            dbg.remove_deadend_nodes();
        });
        eprintln!("[draft/remove_deadends] {}", t);
        eprintln!("[draft] n_nodes={}", dbg.n_nodes());
        eprintln!("[draft] n_edges={}", dbg.n_edges());
        eprintln!("[draft] copy_num_stats_raw={:?}", dbg.copy_num_stats());
        eprintln!("[draft] degree_stats={:?}", dbg.degree_stats());
        // 3
        let freqs = dbg.to_node_freqs() / coverage as f64;
        dbg.set_copy_nums_all_zero();
        let ((copy_nums_approx, cost), t) = timer(|| {
            dbg.min_squared_error_copy_nums_from_freqs_compacted(&freqs, false, &[])
                .unwrap()
        });
        eprintln!("[draft/min_flow] {}", t);
        eprintln!("[draft] approx_cost={}", cost);
        dbg.set_node_copy_nums(&copy_nums_approx);
        eprintln!("[draft] copy_num_stats_approx={:?}", dbg.copy_num_stats());
        dbg
    }
    ///
    ///
    ///
    pub fn create_draft_from_fragment_seqs_with_adjusted_coverage<T>(
        k: usize,
        seqs: T,
        base_coverage: f64,
        ave_read_length: usize,
        p_error: Prob,
        end_node_inference: &EndNodeInference<N::Kmer>,
    ) -> Self
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        let kmer_coverage = kmer_coverage(k, ave_read_length, base_coverage, p_error);
        eprintln!(
            "[draft_frag] k={} ave_read_length={} p_error={} base_coverage={} kmer_coverage={}",
            k, ave_read_length, p_error, base_coverage, kmer_coverage
        );
        Self::create_draft_from_fragment_seqs(k, seqs, kmer_coverage, end_node_inference)
    }
    /// Create draft dbg from fragment reads
    ///
    /// 1. construct dbg without considering starting/ending kmers (= use
    ///    `Dbg::from_fragment_seqs`)
    /// 2. remove 0x and 1x nodes
    /// 3. add starts/ends to all deadend nodes (if Auto) or specified nodes (if Custom)
    /// 4. assign approximate (flow consistent) copy nums by `min_squared_error`.
    ///
    /// ## Known problems
    /// * from fragmented reads, coverage can be overestimated (due to edge effects and k-mer-size
    /// effects)
    ///
    pub fn create_draft_from_fragment_seqs<T>(
        k: usize,
        seqs: T,
        coverage: f64,
        end_node_inference: &EndNodeInference<N::Kmer>,
    ) -> Self
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
        eprintln!("[draft_frag] n_nodes={}", dbg.n_nodes());
        eprintln!("[draft_frag] n_edges={}", dbg.n_edges());
        match end_node_inference {
            EndNodeInference::Auto => {
                dbg.augment_sources_and_sinks();
            }
            EndNodeInference::Custom((starts, ends)) => {
                for (start_kmer, _copy_num) in starts.iter() {
                    dbg.add_starting_kmers(dbg.find_node_from_kmer(start_kmer).unwrap());
                }
                for (end_kmer, _copy_num) in ends.iter() {
                    dbg.add_ending_kmers(dbg.find_node_from_kmer(end_kmer).unwrap());
                }
                dbg.remove_deadend_nodes();
            }
        }
        // 3
        let freqs = dbg.to_node_freqs() / coverage as f64;
        dbg.set_copy_nums_all_zero();
        let fixed_copy_nums = dbg.to_fixed_copy_nums(end_node_inference);
        // debug output
        for (node, copy_num) in fixed_copy_nums.iter() {
            eprintln!("fixed_copy_nums=({}, x{})", node.index(), copy_num);
        }
        let (copy_nums_approx, cost) = dbg
            .min_squared_error_copy_nums_from_freqs_compacted(&freqs, true, &fixed_copy_nums)
            .unwrap();
        eprintln!("[draft] approx_cost={}", cost);
        dbg.set_node_copy_nums(&copy_nums_approx);
        eprintln!("[draft] copy_num_stats_approx={:?}", dbg.copy_num_stats());
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
            |v, weight| {
                if weight.is_emittable() {
                    MinSquaredErrorCopyNumAndFreq::<V2Error>::new(vec![freqs[v]], None, false)
                } else {
                    MinSquaredErrorCopyNumAndFreq::<V2Error>::new(vec![], None, false)
                }
            },
        );
        min_cost_flow_convex_fast(&graph).map(|flow| {
            let cost = total_cost(&graph, &flow);
            // an edge in edbg corresponds to a node in dbg
            // so edgevec for edbg can be converted to nodevec for dbg.
            (flow_to_edgevec(flow).switch_index(), cost)
        })
    }
    ///
    /// by solving min-cost-flow
    /// uses compacted edbg
    ///
    pub fn min_squared_error_copy_nums_from_freqs_compacted(
        &self,
        freqs: &NodeFreqs,
        ignore_startings_and_endings: bool,
        fixed_copy_nums: &[(NodeIndex, CopyNum)],
    ) -> Option<(NodeCopyNums, Cost)> {
        let graph = self.to_compact_edbg_graph();
        // println!("{}", petgraph::dot::Dot::with_config(&graph, &[]));
        let flow_network = graph.map(
            |_, _| (),
            |_, weight| {
                let freqs = weight
                    .origin_edges()
                    .iter()
                    .filter_map(|&edge| {
                        let node = ni(edge.index());
                        let is_target = if ignore_startings_and_endings {
                            self.is_emittable(node) && !self.is_starting_or_ending(node)
                        } else {
                            self.is_emittable(node)
                        };
                        if is_target {
                            Some(freqs[node])
                        } else {
                            None
                        }
                    })
                    .collect();
                // TODO assert that all copy nums specified is consistent?
                let fixed_copy_num = weight.origin_edges().iter().find_map(|edge| {
                    // original node have fixed copy_num?
                    let origin_node = NodeIndex::new(edge.index());
                    fixed_copy_nums
                        .iter()
                        .find(|(node, _)| *node == origin_node)
                        .map(|(_, copy_num)| *copy_num)
                });
                MinSquaredErrorCopyNumAndFreq::<V2Error>::new(freqs, fixed_copy_num, false)
            },
        );
        min_cost_flow_convex_fast(&flow_network).map(|flow| {
            let cost = total_cost(&flow_network, &flow);
            // an edge in edbg corresponds to a node in dbg
            // so edgevec for edbg can be converted to nodevec for dbg.
            // convert flow in compacted-edbg into flow in original edbg
            let flow_in_original = compacted_flow_into_original_flow(self.n_nodes(), &graph, &flow);
            (flow_to_edgevec(flow_in_original).switch_index(), cost)
        })
    }
    ///
    /// get current copy nums as node freqs
    ///
    pub fn to_node_freqs(&self) -> NodeFreqs {
        let copy_nums = self.to_node_copy_nums();
        NodeFreqs::dense_from_vec(
            copy_nums
                .to_inner_vec()
                .into_iter()
                .map(|copy_num| copy_num as Freq)
                .collect(),
            0.0,
        )
    }
    ///
    /// generate `fixed_copy_nums (Vec<(NodeIndex, CopyNum)>)` from `EndNodeInference`.
    ///
    pub fn to_fixed_copy_nums(
        &self,
        end_nodes: &EndNodeInference<N::Kmer>,
    ) -> Vec<(NodeIndex, CopyNum)> {
        match end_nodes {
            EndNodeInference::Auto => Vec::new(),
            EndNodeInference::Custom((starts, ends)) => {
                let mut copy_nums: HashMap<NodeIndex, CopyNum> = HashMap::default();
                for (start, copy_num) in starts {
                    let node = self
                        .find_node_from_kmer(start)
                        .expect("start kmer specified not found!");
                    *copy_nums.entry(node).or_insert(0) += copy_num;
                }
                for (end, copy_num) in ends {
                    let node = self
                        .find_node_from_kmer(end)
                        .expect("end kmer specified not found!");
                    *copy_nums.entry(node).or_insert(0) += copy_num;
                }
                copy_nums.into_iter().collect()
            }
        }
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
    use crate::e2e::{
        generate_difficult_diploid_tandem_repeat_dataset, generate_simple_genome_fragment_dataset,
        generate_simple_genome_mock, generate_tandem_repeat_fragment_dataset,
    };
    use crate::genome;
    use crate::io::cytoscape::NodeAttrVec;
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
        // assert_eq!(cost, 0.0);
        // case 2: all zero freq
        let nf0 = NodeFreqs::new_dense(dbg.n_nodes(), 0.0);
        let (fitted, cost) = dbg.min_squared_error_copy_nums_from_freqs(&nf0).unwrap();
        assert_eq!(fitted.sum(), 0);
        // assert_eq!(cost, 0.0);
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
        // assert!(4.7 <= cost && cost <= 4.8);

        // (2) compact
        let (approx2, cost2) = dbg_raw
            .min_squared_error_copy_nums_from_freqs_compacted(&freq, false, &[])
            .unwrap();
        assert_eq!(approx2.dist(&copy_nums_true), 0);
        assert_eq!(approx2, approx);
        // assert!(4.7 <= cost2 && cost2 <= 4.8);
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
        for (label, dataset) in &[
            ("simple", generate_simple_genome_fragment_dataset()),
            ("tandem_repeat", generate_tandem_repeat_fragment_dataset()),
        ] {
            dataset.show_genome();
            dataset.show_reads();
            println!("coverage={}", dataset.coverage());
            // SimpleDbg::create_draft_from_fragment_seqs(32, dataset.reads(), dataset.coverage());
            let dbg: SimpleDbg<VecKmer> =
                SimpleDbg::create_draft_from_fragment_seqs_with_adjusted_coverage(
                    32,
                    dataset.reads(),
                    dataset.coverage(),
                    50,                       // 50bp read
                    Prob::from_value(0.00_1), // 0.1% error
                    &EndNodeInference::Auto,
                );
            // to check with cytoscape
            let check_with_cytoscape = false;
            if check_with_cytoscape {
                let json = dbg.to_cytoscape();
                write_string(&format!("draft_from_fragment_{}.json", label), &json).unwrap();
            }
            // check if all the true kmers are in the graph
            assert!(dbg.to_copy_nums_of_styled_seqs(dataset.genome()).is_ok());
        }
    }
    #[test]
    fn dbg_create_draft_fragment_bad_approximation_inspection() {
        // "-c 20 -l 100 -p 0.01 --k-init 12 --k-final 100 -U 50 -N 20 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 10 --use-fragment-read";;
        let dataset = generate_difficult_diploid_tandem_repeat_dataset();
        let mut dbg: SimpleDbg<VecKmer> =
            SimpleDbg::create_draft_from_fragment_seqs_with_adjusted_coverage(
                12,
                dataset.reads(),
                dataset.coverage(),
                dataset.reads().average_length(),
                dataset.params().p_error(),
                &EndNodeInference::Auto,
            );
        let (copy_nums_true, _) = dbg.to_copy_nums_of_styled_seqs(dataset.genome()).unwrap();
        let copy_nums_draft = dbg.to_node_copy_nums();
        println!("dist={}", copy_nums_draft.dist(&copy_nums_true));
        println!(
            "max_abs_diff={}",
            copy_nums_draft.max_abs_diff(&copy_nums_true)
        );
        println!("{}", copy_nums_true);
        println!("{}", copy_nums_draft);
        println!("{:?}", copy_nums_draft.diff_element_counts(&copy_nums_true));
        let check_with_cytoscape = false;
        dbg.set_node_copy_nums(&copy_nums_true);
        if check_with_cytoscape {
            let json = dbg.to_cytoscape_with_attrs(&[NodeAttrVec::CopyNum(copy_nums_draft)], &[]);
            write_string(
                &format!("draft_from_fragment_diploid_tandem_repeat.json"),
                &json,
            )
            .unwrap();
        }
    }
}
