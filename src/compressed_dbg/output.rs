//!
//! Output functions for CompressedDBG
//!
use super::CompressedDBG;
use crate::cycles::CycleDirection;
use crate::graph::{Edge, Node, Pos};
use crate::io::cytoscape::Element;
use crate::kmer::kmer::{Kmer, KmerLike};
use crate::prob::Prob;
use crate::stats;
use histo::Histogram;
use std::fmt::Write as FmtWrite;

impl CompressedDBG {
    /// Graphviz dot format
    #[allow(unused_must_use)]
    pub fn as_dot(&self) -> String {
        let mut s = String::new();
        writeln!(&mut s, "digraph cdbg {{");
        for v in self.iter_nodes() {
            // for node
            let kmer = self.kmer(&v);
            writeln!(&mut s, "\t{} [label=\"{}\"];", v.0, kmer);
            // for edges
            for w in self.childs(&v).iter() {
                // writeln!(&mut s, "\t{} -> {} [label=\"{}\"];", v.0, w.0);
                writeln!(&mut s, "\t{} -> {};", v.0, w.0);
            }
        }
        writeln!(&mut s, "}}");
        s
    }
    /// Graphviz dot format
    #[allow(unused_must_use)]
    pub fn as_dot_with_cycle(&self, cycle_id: usize) -> String {
        let cycle = self.cycle_components(cycle_id);
        let mut s = String::new();
        writeln!(&mut s, "digraph cdbg {{");
        for v in self.iter_nodes() {
            // for node
            let kmer = self.kmer(&v);
            if cycle.iter().any(|&(w, _)| w == v) {
                writeln!(&mut s, "\t{} [label=\"{}\" color=red];", v.0, kmer);
            } else {
                writeln!(&mut s, "\t{} [label=\"{}\"];", v.0, kmer);
            }
            // for edges
            for w in self.childs(&v).iter() {
                // writeln!(&mut s, "\t{} -> {} [label=\"{}\"];", v.0, w.0);
                writeln!(&mut s, "\t{} -> {};", v.0, w.0);
            }
        }
        writeln!(&mut s, "}}");
        s
    }
    /// dot with copy_number info on nodes
    #[allow(unused_must_use)]
    pub fn as_dot_with_copy_nums(&self, copy_nums: &[u32]) -> String {
        let mut s = String::new();
        writeln!(&mut s, "digraph cdbg {{");
        for v in self.iter_nodes() {
            // for node
            let kmer = self.kmer(&v);
            let copy_num = copy_nums[v.0];
            writeln!(&mut s, "\t{} [label=\"{} x{}\"];", v.0, kmer, copy_num);
            // for edges
            for w in self.childs(&v).iter() {
                writeln!(&mut s, "\t{} -> {};", v.0, w.0);
            }
        }
        writeln!(&mut s, "}}");
        s
    }
    /// dot with probability (score) on each node
    #[allow(unused_must_use)]
    pub fn as_dot_with_probs(&self, probs: &[Prob]) -> String {
        let mut s = String::new();
        writeln!(&mut s, "digraph cdbg {{");
        for v in self.iter_nodes() {
            // for node
            let kmer = self.kmer(&v);
            let prob = probs[v.0];
            writeln!(&mut s, "\t{} [label=\"{} {}\"];", v.0, kmer, prob);
            // for edges
            for w in self.childs(&v).iter() {
                writeln!(&mut s, "\t{} -> {};", v.0, w.0);
            }
        }
        writeln!(&mut s, "}}");
        s
    }
    /// show a histogram of cycle length distribution
    #[allow(unused_must_use)]
    pub fn as_cycle_histogram(&self) -> String {
        let mut s = String::new();
        writeln!(&mut s, "#cycles={}", self.n_cycles());
        let mut histogram = Histogram::with_buckets(10);
        for i in 0..self.n_cycles() {
            histogram.add(self.cycle_components(i).len() as u64);
        }
        writeln!(&mut s, "{}", histogram);
        s
    }
    // stats related
    pub fn as_dbg_stats(&self) -> stats::DbgStats {
        let mut r = stats::DbgStats::default();
        r.k = self.k;
        r.n_kmers = self.n_kmers as u32;
        r.n_starting_kmers = self.starting_kmers().len() as u32;
        r
    }
    pub fn as_degree_stats(&self) -> stats::DegreeStats {
        let mut r = stats::DegreeStats::default();
        for v in self.iter_nodes() {
            let in_deg = self.parents(&v).len();
            let out_deg = self.childs(&v).len();
            r.in_degs[in_deg] += 1;
            r.out_degs[out_deg] += 1;
        }
        r
    }
    pub fn as_copy_num_stats(&self, copy_nums: &[u32]) -> stats::CopyNumStats {
        let mut r = stats::CopyNumStats::default();
        r.is_consistent = self.is_consistent_copy_num(copy_nums);
        r.total_emitable = self.total_emitable_copy_num(copy_nums);
        r.min = *copy_nums.iter().min().unwrap();
        r.max = *copy_nums.iter().max().unwrap();
        r.total = copy_nums.iter().sum();
        r.average = r.total as f32 / copy_nums.len() as f32;
        r.n_zero_copy_kmer = copy_nums.iter().filter(|c| **c == 0).count() as u32;
        r.n_nonzero_copy_kmer = copy_nums.len() as u32 - r.n_zero_copy_kmer;
        r
    }
    pub fn as_cycle_summary_stats(&self) -> stats::CycleSummaryStats {
        let mut r = stats::CycleSummaryStats::default();
        r.n_cycles = self.n_cycles() as u32;
        let cycle_lengths: Vec<u32> = (0..self.n_cycles())
            .map(|i| self.cycle_components(i).len() as u32)
            .collect();
        if self.n_cycles() == 0 {
            r.min_cycle_len = 0;
            r.max_cycle_len = 0;
            r.average_cycle_len = 0f32;
        } else {
            r.min_cycle_len = *cycle_lengths.iter().min().unwrap();
            r.max_cycle_len = *cycle_lengths.iter().max().unwrap();
            let total: u32 = cycle_lengths.iter().sum();
            r.average_cycle_len = total as f32 / self.n_cycles() as f32;
        }
        r
    }
    pub fn as_cycle_stats(&self, cycle_id: usize) -> stats::CycleStats {
        let mut r = stats::CycleStats::default();
        r.id = cycle_id;
        r.len = self.cycle_components(cycle_id).len();
        r.n_reverse = self
            .cycle_components(cycle_id)
            .iter()
            .map(|(_, dir)| match dir {
                CycleDirection::Reverse => 1,
                _ => 0,
            })
            .sum();
        r
    }
    pub fn as_all_stats(&self, copy_nums: &[u32]) -> stats::AllStats {
        let mut r = stats::AllStats::default();
        r.dbg = Some(self.as_dbg_stats());
        r.copy_num = Some(self.as_copy_num_stats(copy_nums));
        r.degree = Some(self.as_degree_stats());
        r.cycle_summary = Some(self.as_cycle_summary_stats());
        let cycles = (0..self.n_cycles())
            .map(|i| self.as_cycle_stats(i))
            .collect();
        r.cycles = Some(cycles);
        r
    }
    pub fn as_true_kmer_stats(
        &self,
        copy_nums: &[u32],
        copy_nums_true: &[u32],
    ) -> stats::TrueKmerStats {
        let mut r = stats::TrueKmerStats::default();
        r.size = self.total_emitable_copy_num(&copy_nums);
        r.true_size = self.total_emitable_copy_num(&copy_nums_true);
        let true_kmers: Vec<(u32, u32)> = copy_nums
            .iter()
            .zip(copy_nums_true.iter())
            .filter(|(&cn, &cnt)| cn > 0 && cnt > 0)
            .map(|(&cn, &cnt)| (cn, cnt))
            .collect();
        let false_kmers: Vec<(u32, u32)> = copy_nums
            .iter()
            .zip(copy_nums_true.iter())
            .filter(|(&cn, &cnt)| cn > 0 && cnt == 0)
            .map(|(&cn, &cnt)| (cn, cnt))
            .collect();
        r.n_true_kmer = true_kmers.iter().count();
        r.n_false_kmer = false_kmers.iter().count();

        r.true_kmer_max_copy_num = *true_kmers.iter().map(|(cn, _)| cn).max().unwrap_or(&0);
        r.true_kmer_min_copy_num = *true_kmers.iter().map(|(cn, _)| cn).min().unwrap_or(&0);
        r.true_kmer_ave_copy_num =
            true_kmers.iter().map(|(cn, _)| cn).sum::<u32>() as f32 / r.n_true_kmer as f32;

        r.false_kmer_max_copy_num = *false_kmers.iter().map(|(cn, _)| cn).max().unwrap_or(&0);
        r.false_kmer_min_copy_num = *false_kmers.iter().map(|(cn, _)| cn).min().unwrap_or(&0);
        r.false_kmer_ave_copy_num =
            false_kmers.iter().map(|(cn, _)| cn).sum::<u32>() as f32 / r.n_false_kmer as f32;

        r.copy_nums = copy_nums.to_vec();
        r.copy_nums_true = copy_nums_true.to_vec();

        r
    }
    pub fn layout_2d(&self) -> Vec<(Kmer, Pos)> {
        let idg = self.to_indexed_digraph();
        let positions = idg.layout_by_force_atlas2();
        let is_used: Vec<bool> = vec![false; idg.n_nodes()];
        let mut layout = Vec::new();
        for v in self.iter_nodes() {
            let kmer = self.kmer(&v);
            let (v, w) = idg.node_pair(&Edge(v.0));
            // add v
            if !is_used[v.0] {
                layout.push((kmer.prefix(), positions[v.0]))
            }
            // add w
            if !is_used[w.0] {
                layout.push((kmer.suffix(), positions[w.0]))
            }
        }
        layout
    }
    pub fn dump_layout_2d(&self) {
        for (km1mer, pos) in self.layout_2d().iter() {
            println!("N\t{}\t{}\t{}", km1mer, pos.0, pos.1);
        }
    }
    pub fn dump_edge_list(&self, copy_nums: &[u32]) {
        for v in self.iter_nodes() {
            println!("E\t{}\t{}", self.kmer(&v), copy_nums[v.0]);
        }
    }
    /// Create json to feed cytoscape.js graph
    /// format is array of either:
    ///   { group: 'nodes', data: { id } }
    /// or
    ///   { group: 'edges', data: { id, source, target, widths } }
    pub fn to_cytoscape_elements(
        &self,
        copy_nums_list: &[Vec<u32>],
        true_copy_nums: Option<&[u32]>,
    ) -> Vec<Element> {
        let (nodes, edges) = self.to_edge_centric_graph();
        let mut elements = Vec::new();

        // nodes
        let n_nodes = nodes.len();
        for (kmer, node) in nodes.into_iter() {
            elements.push(Element::Node {
                id: node.0,
                label: kmer,
            });
        }

        // edges
        for (i, &(v, w)) in edges.iter().enumerate() {
            let true_width = match true_copy_nums {
                Some(c) => Some(c[i]),
                None => None,
            };
            elements.push(Element::Edge {
                id: i + n_nodes, // element id should be unique, over both nodes and edges
                source: v.0,
                target: w.0,
                label: self.kmer(&Node(i)).clone(),
                widths: copy_nums_list
                    .iter()
                    .map(|copy_nums| copy_nums[i])
                    .collect(),
                true_width,
            });
        }

        elements
    }
    pub fn to_cytoscape_json_with_true(
        &self,
        copy_nums_list: &[Vec<u32>],
        true_copy_nums: Option<&[u32]>,
    ) -> String {
        let elements = self.to_cytoscape_elements(copy_nums_list, true_copy_nums);
        // serde_json::to_string_pretty(&elements).unwrap()
        serde_json::to_string(&elements).unwrap()
    }
    pub fn to_cytoscape_json(&self, copy_nums_list: &[Vec<u32>]) -> String {
        self.to_cytoscape_json_with_true(copy_nums_list, None)
    }
    pub fn to_gexf(&self) -> String {
        let s = String::new();
        s
    }
}
