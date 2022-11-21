//!
//! Definition of Augmented-intersection information collection
//!
//! resolve function
//! * FlowIntersection::resolve_unique
//! * FlowIntersection::resolve_min_flow
//!
//! * is_resolved()
//!
//! * has_freqs
//!
pub mod intersection_graph;
use crate::common::{CopyNum, Freq};
use crate::graph::Bipartite;
use crate::hmmv2::trans_table::EdgeFreqs;
use crate::kmer::kmer::KmerLike;
use crate::min_flow::{min_cost_flow_convex_fast, total_cost, Cost};
use intersection_graph::{IntersectionGraph, IntersectionGraphEdge, IntersectionGraphNode};
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

/*
pub type FlowIntersectionV2<K> = Bipartite<K, FlowIntersectionNode, FlowIntersectionEdge>;

impl<K: KmerLike> FlowIntersectionV2<K> {
    pub fn hoge(&self) {
        println!("hoge");
    }
}
*/

/// Node info
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlowIntersectionNode {
    pub index: NodeIndex,
    pub copy_num: CopyNum,
}

impl FlowIntersectionNode {
    /// Constructor
    pub fn new(index: NodeIndex, copy_num: CopyNum) -> Self {
        FlowIntersectionNode { index, copy_num }
    }
}

impl std::fmt::Display for FlowIntersectionNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}(x{})", self.index.index(), self.copy_num)
    }
}

/// Edge between in-node and out-node
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlowIntersectionEdge {
    pub index: EdgeIndex,
    pub freq: Option<Freq>,
    pub copy_num: Option<CopyNum>,
}

impl FlowIntersectionEdge {
    /// Constructor
    pub fn new(index: EdgeIndex, freq: Option<Freq>, copy_num: Option<CopyNum>) -> Self {
        FlowIntersectionEdge {
            index,
            freq,
            copy_num,
        }
    }
}

impl std::fmt::Display for FlowIntersectionEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.index.index())?;
        match self.freq {
            Some(freq) => write!(f, "(f{})", freq)?,
            None => write!(f, "(f?)")?,
        };
        match self.copy_num {
            Some(copy_num) => write!(f, "(x{})", copy_num),
            None => write!(f, "(x?)"),
        }
    }
}

///
/// Intersection information corresponding to a k-1-mer.
///
/// in node-centric de bruijn graph, a k-1-mer can have
/// incoming kmers and outgoing kmers (less than 4 nodes).
///
/// km1mer:    XXX
/// in_nodes:  AXXX, CXXX, GXXX, TXXX
/// out_nodes: XXXA, XXXC, XXXG, XXXT
///
#[derive(Clone, Debug)]
pub struct FlowIntersection<K: KmerLike> {
    pub bi: Bipartite<K, FlowIntersectionNode, FlowIntersectionEdge>,
}

impl<K: KmerLike> FlowIntersection<K> {
    pub fn new(
        km1mer: K,
        in_nodes: Vec<FlowIntersectionNode>,
        out_nodes: Vec<FlowIntersectionNode>,
        edges: Vec<FlowIntersectionEdge>,
    ) -> Self {
        FlowIntersection {
            bi: Bipartite::new(km1mer, in_nodes, out_nodes, edges),
        }
    }
    pub fn from<F: Fn(usize, usize) -> FlowIntersectionEdge>(
        km1mer: K,
        in_nodes: Vec<FlowIntersectionNode>,
        out_nodes: Vec<FlowIntersectionNode>,
        edge_fn: F,
    ) -> Self {
        FlowIntersection {
            bi: Bipartite::from(km1mer, in_nodes, out_nodes, edge_fn),
        }
    }
}

///
/// Property getters
///
impl<K: KmerLike> FlowIntersection<K> {
    pub fn n_in_nodes(&self) -> usize {
        self.bi.n_in()
    }
    pub fn n_out_nodes(&self) -> usize {
        self.bi.n_out()
    }
    pub fn km1mer(&self) -> &K {
        &self.bi.id
    }
    pub fn iter_in_node_indexes(&self) -> impl Iterator<Item = NodeIndex> + '_ {
        self.bi.in_nodes.iter().map(move |v| v.index)
    }
    pub fn iter_out_node_indexes(&self) -> impl Iterator<Item = NodeIndex> + '_ {
        self.bi.out_nodes.iter().map(move |v| v.index)
    }
    pub fn iter_in_nodes(&self) -> impl Iterator<Item = &FlowIntersectionNode> + '_ {
        self.bi.in_nodes.iter()
    }
    pub fn iter_out_nodes(&self) -> impl Iterator<Item = &FlowIntersectionNode> + '_ {
        self.bi.out_nodes.iter()
    }
    pub fn iter_edges(&self) -> impl Iterator<Item = (usize, usize, &FlowIntersectionEdge)> + '_ {
        self.bi.iter_edges()
    }
    pub fn in_node_index(&self, index: usize) -> NodeIndex {
        self.bi.in_node(index).index
    }
    pub fn out_node_index(&self, index: usize) -> NodeIndex {
        self.bi.out_node(index).index
    }
    pub fn in_node(&self, index: usize) -> &FlowIntersectionNode {
        self.bi.in_node(index)
    }
    pub fn out_node(&self, index: usize) -> &FlowIntersectionNode {
        self.bi.out_node(index)
    }
    pub fn edge(&self, i: usize, j: usize) -> &FlowIntersectionEdge {
        self.bi.edge(i, j)
    }
    pub fn in_node_copy_nums(&self) -> Vec<CopyNum> {
        self.bi.in_nodes.iter().map(move |v| v.copy_num).collect()
    }
    pub fn out_node_copy_nums(&self) -> Vec<CopyNum> {
        self.bi.out_nodes.iter().map(move |v| v.copy_num).collect()
    }
}

///
/// Appending freq information
///
impl<K: KmerLike> FlowIntersection<K> {
    ///
    /// Append freq information to edges using EdgeFreq vector.
    ///
    pub fn augment_freqs(mut self, freqs: &EdgeFreqs) -> Self {
        for edge in self.bi.edges.iter_mut() {
            edge.freq = Some(freqs[edge.index]);
        }
        self
    }
    ///
    /// Check if each edge has freq information or not.
    ///
    pub fn has_freqs(&self) -> bool {
        self.iter_edges().all(|(_, _, e)| e.freq.is_some())
    }
    ///
    /// sum of freqs of all edges between in/out nodes
    ///
    pub fn total_freq(&self) -> Option<Freq> {
        if self.has_freqs() {
            let freq: Freq = self.iter_edges().map(|(_, _, e)| e.freq.unwrap()).sum();
            Some(freq)
        } else {
            None
        }
    }
}

///
/// Upconvert related
///
impl<K: KmerLike> FlowIntersection<K> {
    ///
    /// sum of copy nums of in_nodes
    ///
    pub fn total_in_copy_num(&self) -> CopyNum {
        self.iter_in_nodes().map(|n| n.copy_num).sum()
    }
    ///
    /// sum of copy nums of out_nodes
    ///
    pub fn total_out_copy_num(&self) -> CopyNum {
        self.iter_out_nodes().map(|n| n.copy_num).sum()
    }
    ///
    /// ambiguous <=> not uniquely resolvable and not tip_intersection
    ///
    pub fn is_ambiguous(&self) -> bool {
        !self.is_tip_intersection() && !self.can_unique_resolvable()
    }
    ///
    /// sum of copynums of in_nodes/out_nodes are the same?
    ///
    pub fn has_valid_node_copy_nums(&self) -> bool {
        self.total_in_copy_num() == self.total_out_copy_num()
    }
    pub fn can_unique_resolvable(&self) -> bool {
        self.n_in_nodes() == 1 || self.n_out_nodes() == 1
    }
    pub fn all_edges_has_copy_num(&self) -> bool {
        self.iter_edges().all(|(_, _, e)| e.copy_num.is_some())
    }
    pub fn all_edge_copy_nums_consistent(&self) -> bool {
        let n_in = self.n_in_nodes();
        let n_out = self.n_out_nodes();
        // for all in_node, check that the sum-copy-num from the in_node should be the same as
        // copy-num of the node
        for in_node in 0..n_in {
            let sum_edge_copy_nums: CopyNum = (0..n_out)
                .map(|out_node| self.edge(in_node, out_node).copy_num.unwrap())
                .sum();
            let in_node_copy_num = self.in_node(in_node).copy_num;
            if sum_edge_copy_nums != in_node_copy_num {
                return false;
            }
        }
        // also for out_node
        for out_node in 0..n_out {
            let sum_edge_copy_nums: CopyNum = (0..n_in)
                .map(|in_node| self.edge(in_node, out_node).copy_num.unwrap())
                .sum();
            let out_node_copy_num = self.out_node(out_node).copy_num;
            if sum_edge_copy_nums != out_node_copy_num {
                return false;
            }
        }
        // if it passed all checks, the copy nums are consistent.
        true
    }
    ///
    /// Check if the km1mer of the intersection is NNNN.
    /// If it is a tip intersection, the edge copy number should not be resolved.
    ///
    pub fn is_tip_intersection(&self) -> bool {
        self.km1mer().is_null()
    }
    ///
    /// Check if the km1mer of the intersection is ending (e.g. XNNN or NNYY).
    /// Ending intersection should be a simple (n_in==1 or n_out==1).
    ///
    /// ```text
    /// YXNNN -->
    ///           XNNNN
    /// ZXNNN -->
    /// ```
    ///
    /// node(`YXNNN`) should have only one child(`XNNNN`), so the intersection `XNNN` should be
    /// simple.
    ///
    pub fn is_end_intersection(&self) -> bool {
        let km1mer = self.km1mer();
        !km1mer.is_null() && km1mer.has_null()
    }
    ///
    /// FlowIntersection is resolved or not.
    ///
    /// * For a normal intersection, it is resolved if the all edges has its own copy number.
    /// * For a tip intersection, the edge copy number cannot be assigned, so it will always be marked
    /// as resolved.
    ///
    pub fn is_resolved(&self) -> bool {
        if self.is_tip_intersection() {
            true
        } else {
            self.all_edges_has_copy_num()
        }
    }
    /// Do appropriate conversion.
    ///
    /// * if uniquly convertable, do a simple conversion
    /// * otherwise, do a optimize conversion using min flow.
    ///
    /// Score of the min-flow is returned when min_flow is used.
    pub fn resolve(&self) -> (FlowIntersection<K>, Option<Cost>) {
        if self.can_unique_resolvable() {
            (self.resolve_unique(), None)
        } else {
            let (fio, cost) = self.resolve_min_flow();
            (fio, Some(cost))
        }
    }
    /// Get copy-number-resolved FlowIntersection by applying unique (obvious)
    /// conversion.
    ///
    /// If `n_in == 1 or n_out == 1`, a copy number of each edge `in -> out`
    /// should have the same copy number of in/out.
    ///
    pub fn resolve_unique(&self) -> FlowIntersection<K> {
        assert!(self.can_unique_resolvable());
        let mut opt = self.clone();

        if self.n_in_nodes() == 1 {
            for i in 0..self.n_out_nodes() {
                opt.bi.edges[i].copy_num = Some(opt.bi.out_node(i).copy_num);
            }
        } else if self.n_out_nodes() == 1 {
            for i in 0..self.n_in_nodes() {
                opt.bi.edges[i].copy_num = Some(opt.bi.in_node(i).copy_num);
            }
        }

        // check if all edges have its copy num
        assert!(opt.all_edges_has_copy_num());
        opt
    }
    ///
    ///
    ///
    pub fn resolve_naive(&self) -> FlowIntersection<K> {
        let mut opt = self.clone();
        let n_in = opt.n_in_nodes();
        let n_out = opt.n_out_nodes();
        // initialize with zero copy num
        for in_node in 0..n_in {
            for out_node in 0..n_out {
                opt.bi.edge_mut(in_node, out_node).copy_num = Some(0);
            }
        }
        let mut in_node_copy_nums = opt.in_node_copy_nums();
        let mut out_node_copy_nums = opt.out_node_copy_nums();
        // check if there is remaining copy_nums in node
        let is_remaining = |in_node_copy_nums: &[usize], out_node_copy_nums: &[usize]| {
            let in_remaining: usize = in_node_copy_nums.iter().sum();
            let out_remaining: usize = out_node_copy_nums.iter().sum();
            in_remaining > 0 || out_remaining > 0
        };
        while is_remaining(&in_node_copy_nums, &out_node_copy_nums) {
            for in_node in 0..n_in {
                for out_node in 0..n_out {
                    if in_node_copy_nums[in_node] > 0 && out_node_copy_nums[out_node] > 0 {
                        // add one to this edge
                        let c = opt.bi.edge(in_node, out_node).copy_num.unwrap();
                        opt.bi.edge_mut(in_node, out_node).copy_num = Some(c + 1);
                        in_node_copy_nums[in_node] -= 1;
                        out_node_copy_nums[out_node] -= 1;
                    }
                }
            }
        }
        // TODO check that all copy nums on edges are consistent
        assert!(opt.all_edges_has_copy_num());
        assert!(opt.all_edge_copy_nums_consistent());
        opt
    }
    ///
    /// Get optimized copy numbers of edges and its min-flow cost.
    /// by converting the bipartite into flow network definitions
    ///
    pub fn resolve_min_flow(&self) -> (FlowIntersection<K>, Cost) {
        let flow_graph = self.to_flow_graph();
        match min_cost_flow_convex_fast(&flow_graph) {
            Some(flow) => {
                let mut opt = self.clone();

                for e in flow_graph.edge_indices() {
                    let ew = flow_graph.edge_weight(e).unwrap();
                    match ew {
                        IntersectionGraphEdge::Spanning {
                            freq,
                            edge_index,
                            edge_index_in_intersection,
                        } => {
                            // modify copy_num in FlowIntersectionEdge by optimized one
                            let copy_num = flow[e];
                            opt.bi.edges[*edge_index_in_intersection].copy_num = Some(copy_num);
                        }
                        _ => {}
                    };
                }
                let cost = total_cost(&flow_graph, &flow);
                assert!(opt.all_edges_has_copy_num());
                (opt, cost)
            }
            None => {
                panic!("flow not found");
            }
        }
    }
    ///
    /// Convert FlowIntersection into DiGraph
    /// that can be solved as a convex min flow problem.
    ///
    pub fn to_flow_graph(&self) -> IntersectionGraph {
        let mut g = IntersectionGraph::new();

        // (1) Node
        // (1-a) source and terminal
        let s = g.add_node(IntersectionGraphNode::Source);
        let t = g.add_node(IntersectionGraphNode::Terminal);
        // (1-b) left/right nodes
        let vs: Vec<NodeIndex> = self
            .bi
            .iter_in_nodes()
            .map(|v| g.add_node(IntersectionGraphNode::Left(v.index)))
            .collect();
        let ws: Vec<NodeIndex> = self
            .bi
            .iter_out_nodes()
            .map(|w| g.add_node(IntersectionGraphNode::Right(w.index)))
            .collect();

        // (2) Edge
        // (2-a) LoopBack
        g.add_edge(t, s, IntersectionGraphEdge::LoopBack);
        // (2-b) Node
        for (i, n) in self.bi.iter_in_nodes().enumerate() {
            g.add_edge(
                s,
                vs[i],
                IntersectionGraphEdge::Node {
                    copy_num: n.copy_num,
                },
            );
        }
        for (i, n) in self.bi.iter_out_nodes().enumerate() {
            g.add_edge(
                ws[i],
                t,
                IntersectionGraphEdge::Node {
                    copy_num: n.copy_num,
                },
            );
        }
        // (2-c) Spanning
        for i in 0..self.bi.n_in() {
            for j in 0..self.bi.n_out() {
                let e = self.bi.edge(i, j);
                let ei = self.bi.edge_index(i, j);
                g.add_edge(
                    vs[i],
                    ws[j],
                    IntersectionGraphEdge::Spanning {
                        freq: e
                            .freq
                            .expect("intersection without freqs cannot be converted to flow graph"),
                        edge_index: e.index,
                        edge_index_in_intersection: ei,
                    },
                );
            }
        }

        g
    }
    /// (For debug)
    /// convert into copy numbers list
    ///
    /// TODO convert `Vec<Option<T>>` into `Option<Vec<T>>`.
    ///
    fn to_edge_copy_nums(&self) -> Vec<Option<CopyNum>> {
        self.iter_edges().map(|(_, _, e)| e.copy_num).collect()
    }
}

impl<K: KmerLike> std::fmt::Display for FlowIntersection<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "FlowIntersection\n{}", self.bi)
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};
    use crate::kmer::veckmer::VecKmer;
    use crate::min_flow::utils::clamped_log;
    use petgraph::dot::Dot;

    #[test]
    fn flow_intersection_construction() {
        let in_nodes = vec![
            FlowIntersectionNode::new(ni(10), 2),
            FlowIntersectionNode::new(ni(11), 4),
        ];
        let out_nodes = vec![
            FlowIntersectionNode::new(ni(13), 1),
            FlowIntersectionNode::new(ni(14), 5),
        ];
        let edges = vec![
            FlowIntersectionEdge::new(ei(20), Some(5.1), None),
            FlowIntersectionEdge::new(ei(21), Some(5.0), None),
            FlowIntersectionEdge::new(ei(22), Some(4.9), None),
            FlowIntersectionEdge::new(ei(23), Some(4.8), None),
        ];
        let edge_copy_nums_true = vec![1, 1, 0, 4];
        let fi = FlowIntersection::new(VecKmer::from_bases(b"TCG"), in_nodes, out_nodes, edges);
        println!("{}", fi);
        let g = fi.to_flow_graph();
        println!("{:?}", Dot::with_config(&g, &[]));

        let (fi_opt, cost) = fi.resolve();
        println!("{}", fi_opt);
        println!("cost={:?}", cost);
        let cost_true = 5.1 * clamped_log(1)
            + 5.0 * clamped_log(1)
            + 4.9 * clamped_log(0)
            + 4.8 * clamped_log(4);
        println!("cost_true={}", cost_true);
        assert_eq!(cost, Some(-cost_true));
        for (i, e) in fi_opt.bi.edges.iter().enumerate() {
            assert_eq!(e.index, fi.bi.edges[i].index);
            assert_eq!(e.freq, fi.bi.edges[i].freq);
            assert_eq!(e.copy_num, Some(edge_copy_nums_true[i]));
        }
    }
    #[test]
    fn flow_intersection_one() {
        let in_nodes = vec![
            FlowIntersectionNode::new(ni(10), 1),
            FlowIntersectionNode::new(ni(11), 1),
        ];
        let out_nodes = vec![
            FlowIntersectionNode::new(ni(13), 1),
            FlowIntersectionNode::new(ni(14), 1),
        ];
        let edges = vec![
            FlowIntersectionEdge::new(ei(20), Some(0.9), None),
            FlowIntersectionEdge::new(ei(21), Some(0.0), None),
            FlowIntersectionEdge::new(ei(22), Some(0.0), None),
            FlowIntersectionEdge::new(ei(23), Some(1.1), None),
        ];
        let kmer = VecKmer::from_bases(b"TCG");
        let fi = FlowIntersection::new(kmer, in_nodes, out_nodes, edges);
        let (fi_opt, cost) = fi.resolve();
        println!("{}", fi_opt);
        assert_eq!(
            fi_opt.to_edge_copy_nums(),
            vec![Some(1), Some(0), Some(0), Some(1)]
        );

        let in_nodes = vec![
            FlowIntersectionNode::new(ni(10), 1),
            FlowIntersectionNode::new(ni(11), 1),
        ];
        let out_nodes = vec![
            FlowIntersectionNode::new(ni(13), 1),
            FlowIntersectionNode::new(ni(14), 1),
        ];
        let edges = vec![
            FlowIntersectionEdge::new(ei(20), Some(0.1), None),
            FlowIntersectionEdge::new(ei(21), Some(1.1), None),
            FlowIntersectionEdge::new(ei(22), Some(0.5), None),
            FlowIntersectionEdge::new(ei(23), Some(0.5), None),
        ];
        let kmer = VecKmer::from_bases(b"TCG");
        let fi = FlowIntersection::new(kmer, in_nodes, out_nodes, edges);
        let (fi_opt, cost) = fi.resolve();
        println!("{}", fi_opt);
        println!("cost={:?}", cost);
        let cost_true = 0.1 * clamped_log(0)
            + 1.1 * clamped_log(1)
            + 0.5 * clamped_log(1)
            + 0.5 * clamped_log(0);
        println!("cost_true={}", cost_true);
        assert_eq!(cost, Some(-cost_true));
        assert_eq!(
            fi_opt.to_edge_copy_nums(),
            vec![Some(0), Some(1), Some(1), Some(0)]
        );
    }
    #[test]
    fn flow_intersection_convert_unique() {
        // (0)
        let in_nodes = vec![FlowIntersectionNode::new(ni(0), 5)];
        let out_nodes = vec![FlowIntersectionNode::new(ni(1), 2)];
        let edges = vec![FlowIntersectionEdge::new(ei(0), Some(0.9), None)];
        let kmer = VecKmer::from_bases(b"TCG");
        let fi = FlowIntersection::new(kmer, in_nodes, out_nodes, edges);
        assert!(!fi.has_valid_node_copy_nums());

        // (1) obviously converable case (n_in = 1)
        let in_nodes = vec![FlowIntersectionNode::new(ni(0), 5)];
        let out_nodes = vec![
            FlowIntersectionNode::new(ni(1), 2),
            FlowIntersectionNode::new(ni(2), 3),
        ];
        let edges = vec![
            FlowIntersectionEdge::new(ei(0), Some(0.9), None),
            FlowIntersectionEdge::new(ei(1), Some(0.0), None),
        ];
        let kmer = VecKmer::from_bases(b"TCG");
        let fi = FlowIntersection::new(kmer, in_nodes, out_nodes, edges);
        println!("{}", fi);
        assert!(!fi.all_edges_has_copy_num());
        assert!(fi.has_valid_node_copy_nums());
        assert!(fi.can_unique_resolvable());
        let fio = fi.resolve_unique();
        println!("{}", fio);
        assert_eq!(
            fio.bi.edges,
            vec![
                FlowIntersectionEdge::new(ei(0), Some(0.9), Some(2)),
                FlowIntersectionEdge::new(ei(1), Some(0.0), Some(3)),
            ]
        );
        let (fio2, cost) = fi.resolve();
        assert_eq!(fio.bi.edges, fio2.bi.edges);

        // (3) obviously converable case (n_out = 1)
        let in_nodes = vec![
            FlowIntersectionNode::new(ni(0), 5),
            FlowIntersectionNode::new(ni(1), 3),
        ];
        let out_nodes = vec![FlowIntersectionNode::new(ni(2), 8)];
        let edges = vec![
            FlowIntersectionEdge::new(ei(0), Some(0.9), None),
            FlowIntersectionEdge::new(ei(1), Some(0.0), None),
        ];
        let kmer = VecKmer::from_bases(b"TCG");
        let fi = FlowIntersection::new(kmer, in_nodes, out_nodes, edges);
        println!("{}", fi);
        assert!(!fi.all_edges_has_copy_num());
        assert!(fi.has_valid_node_copy_nums());
        assert!(fi.can_unique_resolvable());
        let fio = fi.resolve_unique();
        println!("{}", fio);
        assert_eq!(
            fio.bi.edges,
            vec![
                FlowIntersectionEdge::new(ei(0), Some(0.9), Some(5)),
                FlowIntersectionEdge::new(ei(1), Some(0.0), Some(3)),
            ]
        );
        let (fio2, cost) = fi.resolve();
        assert_eq!(fio.bi.edges, fio2.bi.edges);
    }
    #[test]
    fn flow_intersection_convert_naive_01() {
        // [0] multiple node case
        let in_nodes = vec![
            FlowIntersectionNode::new(ni(0), 5),
            FlowIntersectionNode::new(ni(1), 3),
        ];
        let out_nodes = vec![
            FlowIntersectionNode::new(ni(10), 3),
            FlowIntersectionNode::new(ni(11), 2),
            FlowIntersectionNode::new(ni(12), 3),
        ];
        let edges = vec![
            FlowIntersectionEdge::new(ei(0), None, None),
            FlowIntersectionEdge::new(ei(1), None, None),
            FlowIntersectionEdge::new(ei(2), None, None),
            FlowIntersectionEdge::new(ei(3), None, None),
            FlowIntersectionEdge::new(ei(4), None, None),
            FlowIntersectionEdge::new(ei(5), None, None),
        ];
        let kmer = VecKmer::from_bases(b"TCG");
        let fi = FlowIntersection::new(kmer, in_nodes, out_nodes, edges);
        println!("{}", fi);
        let fio = fi.resolve_naive();
        println!("{}", fio);
        assert_eq!(
            fio.to_edge_copy_nums(),
            vec![Some(2), Some(1), Some(2), Some(1), Some(1), Some(1)]
        );
    }
    #[test]
    fn flow_intersection_convert_naive_02() {
        // [0] simple node case
        let in_nodes = vec![
            FlowIntersectionNode::new(ni(0), 1),
            FlowIntersectionNode::new(ni(1), 1),
        ];
        let out_nodes = vec![
            FlowIntersectionNode::new(ni(10), 1),
            FlowIntersectionNode::new(ni(11), 1),
        ];
        let edges = vec![
            FlowIntersectionEdge::new(ei(0), None, None),
            FlowIntersectionEdge::new(ei(1), None, None),
            FlowIntersectionEdge::new(ei(2), None, None),
            FlowIntersectionEdge::new(ei(3), None, None),
        ];
        let kmer = VecKmer::from_bases(b"TCG");
        let fi = FlowIntersection::new(kmer, in_nodes, out_nodes, edges);
        println!("{}", fi);
        let fio = fi.resolve_naive();
        println!("{}", fio);
        assert_eq!(
            fio.to_edge_copy_nums(),
            vec![Some(1), Some(0), Some(0), Some(1)]
        );
    }
}
