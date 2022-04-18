//!
//! Definition of Augmented-intersection information collection
//!
use super::intersection_graph::{IntersectionGraph, IntersectionGraphEdge, IntersectionGraphNode};
use crate::common::{CopyNum, Freq};
use crate::graph::Bipartite;
use crate::kmer::kmer::KmerLike;
use crate::min_flow::{min_cost_flow_convex_fast, total_cost, Cost};
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

pub type FlowIntersectionV2<K> = Bipartite<K, FlowIntersectionNode, FlowIntersectionEdge>;

impl<K: KmerLike> FlowIntersectionV2<K> {
    pub fn hoge(&self) {
        println!("hoge");
    }
}

/// Node info
#[derive(Clone, Debug, PartialEq)]
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
#[derive(Clone, Debug, PartialEq)]
pub struct FlowIntersectionEdge {
    pub index: EdgeIndex,
    pub freq: Freq,
    pub copy_num: Option<CopyNum>,
}

impl FlowIntersectionEdge {
    /// Constructor
    pub fn new(index: EdgeIndex, freq: Freq, copy_num: Option<CopyNum>) -> Self {
        FlowIntersectionEdge {
            index,
            freq,
            copy_num,
        }
    }
}

impl std::fmt::Display for FlowIntersectionEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}({})", self.index.index(), self.freq)?;
        match self.copy_num {
            Some(copy_num) => write!(f, "(x{})", copy_num),
            None => write!(f, "(x?)"),
        }
    }
}

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
}

///
/// Upconvert related
///
impl<K: KmerLike> FlowIntersection<K> {
    ///
    /// sum of copynums of in_nodes/out_nodes are the same?
    ///
    pub fn has_valid_node_copy_nums(&self) -> bool {
        let sum_in_nodes: usize = self.bi.iter_in_nodes().map(|n| n.copy_num).sum();
        let sum_out_nodes: usize = self.bi.iter_out_nodes().map(|n| n.copy_num).sum();
        sum_in_nodes == sum_out_nodes
    }
    pub fn can_uniquely_convertable(&self) -> bool {
        self.n_in_nodes() == 1 || self.n_out_nodes() == 1
    }
    pub fn all_edges_has_copy_num(&self) -> bool {
        self.bi.iter_edges().all(|(_, _, e)| e.copy_num.is_some())
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
    /// * if uniquly convertable, do a simple conversion
    /// * otherwise, do a optimize conversion using min flow.
    ///
    /// Score of the min-flow is returned when min_flow is used.
    pub fn convert(&self) -> (FlowIntersection<K>, Option<Cost>) {
        if self.can_uniquely_convertable() {
            (self.unique_convert(), None)
        } else {
            let (fio, cost) = self.optimize();
            (fio, Some(cost))
        }
    }
    /// Get copy-number-resolved FlowIntersection by applying unique (obvious)
    /// conversion.
    ///
    /// If `n_in == 1 or n_out == 1`, a copy number of each edge `in -> out`
    /// should have the same copy number of in/out.
    ///
    pub fn unique_convert(&self) -> FlowIntersection<K> {
        assert!(self.can_uniquely_convertable());
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
    /// Get optimized copy numbers of edges and its min-flow cost.
    /// by converting the bipartite into flow network definitions
    ///
    pub fn optimize(&self) -> (FlowIntersection<K>, Cost) {
        println!("optimizing: {}", self);
        let flow_graph = self.to_flow_graph();
        match min_cost_flow_convex_fast(&flow_graph) {
            Some(flow) => {
                println!("flow_found = {}", flow);
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
                println!("optimized: {}", opt);
                println!("cost: {}", cost);
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
                        freq: e.freq,
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
        self.bi.edges.iter().map(|e| e.copy_num).collect()
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
            FlowIntersectionEdge::new(ei(20), 5.1, None),
            FlowIntersectionEdge::new(ei(21), 5.0, None),
            FlowIntersectionEdge::new(ei(22), 4.9, None),
            FlowIntersectionEdge::new(ei(23), 4.8, None),
        ];
        let edge_copy_nums_true = vec![1, 1, 0, 4];
        let fi = FlowIntersection::new(VecKmer::from_bases(b"TCG"), in_nodes, out_nodes, edges);
        println!("{}", fi);
        let g = fi.to_flow_graph();
        println!("{:?}", Dot::with_config(&g, &[]));

        let (fi_opt, cost) = fi.optimize();
        println!("{}", fi_opt);
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
            FlowIntersectionEdge::new(ei(20), 0.9, None),
            FlowIntersectionEdge::new(ei(21), 0.0, None),
            FlowIntersectionEdge::new(ei(22), 0.0, None),
            FlowIntersectionEdge::new(ei(23), 1.1, None),
        ];
        let kmer = VecKmer::from_bases(b"TCG");
        let fi = FlowIntersection::new(kmer, in_nodes, out_nodes, edges);
        let (fi_opt, cost) = fi.optimize();
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
            FlowIntersectionEdge::new(ei(20), 0.1, None),
            FlowIntersectionEdge::new(ei(21), 1.1, None),
            FlowIntersectionEdge::new(ei(22), 0.5, None),
            FlowIntersectionEdge::new(ei(23), 0.5, None),
        ];
        let kmer = VecKmer::from_bases(b"TCG");
        let fi = FlowIntersection::new(kmer, in_nodes, out_nodes, edges);
        let (fi_opt, cost) = fi.optimize();
        println!("{}", fi_opt);
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
        let edges = vec![FlowIntersectionEdge::new(ei(0), 0.9, None)];
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
            FlowIntersectionEdge::new(ei(0), 0.9, None),
            FlowIntersectionEdge::new(ei(1), 0.0, None),
        ];
        let kmer = VecKmer::from_bases(b"TCG");
        let fi = FlowIntersection::new(kmer, in_nodes, out_nodes, edges);
        println!("{}", fi);
        assert!(!fi.all_edges_has_copy_num());
        assert!(fi.has_valid_node_copy_nums());
        assert!(fi.can_uniquely_convertable());
        let fio = fi.unique_convert();
        println!("{}", fio);
        assert_eq!(
            fio.bi.edges,
            vec![
                FlowIntersectionEdge::new(ei(0), 0.9, Some(2)),
                FlowIntersectionEdge::new(ei(1), 0.0, Some(3)),
            ]
        );
        let (fio2, cost) = fi.convert();
        assert_eq!(fio.bi.edges, fio2.bi.edges);

        // (3) obviously converable case (n_out = 1)
        let in_nodes = vec![
            FlowIntersectionNode::new(ni(0), 5),
            FlowIntersectionNode::new(ni(1), 3),
        ];
        let out_nodes = vec![FlowIntersectionNode::new(ni(2), 8)];
        let edges = vec![
            FlowIntersectionEdge::new(ei(0), 0.9, None),
            FlowIntersectionEdge::new(ei(1), 0.0, None),
        ];
        let kmer = VecKmer::from_bases(b"TCG");
        let fi = FlowIntersection::new(kmer, in_nodes, out_nodes, edges);
        println!("{}", fi);
        assert!(!fi.all_edges_has_copy_num());
        assert!(fi.has_valid_node_copy_nums());
        assert!(fi.can_uniquely_convertable());
        let fio = fi.unique_convert();
        println!("{}", fio);
        assert_eq!(
            fio.bi.edges,
            vec![
                FlowIntersectionEdge::new(ei(0), 0.9, Some(5)),
                FlowIntersectionEdge::new(ei(1), 0.0, Some(3)),
            ]
        );
        let (fio2, cost) = fi.convert();
        assert_eq!(fio.bi.edges, fio2.bi.edges);
    }
}
