//!
//! `GenomeGraph`
//! Node: linear sequence
//! Edge: its adjacency
//!
use super::seq_graph::{get_start_points, SimpleSeqEdge, SimpleSeqGraph, SimpleSeqNode};
use crate::common::{CopyNum, PositionedReads, PositionedSequence, Reads, Seq, Sequence};
use crate::graph::seq_graph::SeqGraph;
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::sample::{ReadAmount, SampleProfile, StartPoints};
use itertools::Itertools;
use petgraph::dot::Dot;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use petgraph::visit::IntoNodeReferences;
use petgraph::Direction;
use std::collections::HashMap;

/// Read sampling profile
///
#[derive(Clone, Debug)]
pub struct ReadProfile {
    pub sample_profile: SampleProfile,
    pub phmm_params: PHMMParams,
}

/// GenomeGraph
pub struct GenomeGraph(pub DiGraph<GenomeNode, GenomeEdge>);

/// Node in GenomeGraph is a copy-number-assigned sequence fragment
pub struct GenomeNode {
    /// sequence fragment of this node
    pub seq: Sequence,
    /// copy number of the sequence
    pub copy_num: CopyNum,
}

impl GenomeNode {
    pub fn new(seq: &[u8], copy_num: CopyNum) -> Self {
        GenomeNode {
            seq: seq.to_vec(),
            copy_num,
        }
    }
}

pub struct GenomeEdge {
    /// copy number of the transition
    /// it assumes random transition if `None`
    pub copy_num: Option<CopyNum>,
}

impl GenomeEdge {
    pub fn new(copy_num: Option<CopyNum>) -> Self {
        GenomeEdge { copy_num }
    }
}

#[derive(Clone, Debug, Copy, PartialEq)]
pub struct GenomeGraphPos {
    /// index of GenomeNode
    node: NodeIndex,
    /// position in GenomeNode
    pos: usize,
}

impl GenomeGraphPos {
    /// Constructor
    pub fn new(node: NodeIndex, pos: usize) -> Self {
        GenomeGraphPos { node, pos }
    }
    /// is first base or not (pos == 0)
    pub fn is_first_base(&self) -> bool {
        self.pos == 0
    }
}

impl std::fmt::Display for GenomeGraphPos {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({},{})", self.node.index(), self.pos)
    }
}

pub type GenomeGraphPosVec = Vec<GenomeGraphPos>;

impl std::fmt::Display for GenomeNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let seq = std::str::from_utf8(&self.seq).unwrap();
        write!(f, "{} (x{})", seq, self.copy_num)
    }
}

impl std::fmt::Display for GenomeEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.copy_num {
            Some(copy_num) => {
                write!(f, "x{}", copy_num)
            }
            None => {
                write!(f, "")
            }
        }
    }
}

impl std::fmt::Display for GenomeGraph {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.0, &[]))
    }
}

impl GenomeGraph {
    ///
    /// Create `GenomeGraph` from the sequences
    ///
    pub fn from_seqs(seqs: &[Sequence]) -> Self {
        let mut g = DiGraph::new();
        for seq in seqs {
            g.add_node(GenomeNode::new(seq, 1));
        }
        GenomeGraph(g)
    }
    /// Get the number of nodes in the GenomeGraph
    pub fn node_count(&self) -> usize {
        self.0.node_count()
    }
    /// Get the number of edges in the GenomeGraph
    pub fn edge_count(&self) -> usize {
        self.0.edge_count()
    }
    /// wrapper of DiGraph::node_indices
    pub fn node_indices(&self) -> impl Iterator<Item = NodeIndex> {
        self.0.node_indices()
    }
    /// wrapper of DiGraph::edge_indices
    pub fn edge_indices(&self) -> impl Iterator<Item = EdgeIndex> {
        self.0.edge_indices()
    }
    /// wrapper of DiGraph::node_weight
    pub fn node_weight(&self, node: NodeIndex) -> Option<&GenomeNode> {
        self.0.node_weight(node)
    }
    /// wrapper of DiGraph::edge_weight
    pub fn edge_weight(&self, edge: EdgeIndex) -> Option<&GenomeEdge> {
        self.0.edge_weight(edge)
    }
    /// wrapper of DiGraph::edge_endpoints
    pub fn edge_endpoints(&self, edge: EdgeIndex) -> Option<(NodeIndex, NodeIndex)> {
        self.0.edge_endpoints(edge)
    }
    /// if node has no incoming edge, it is a start point node
    pub fn is_start_point_node(&self, node: NodeIndex) -> bool {
        self.0.edges_directed(node, Direction::Incoming).count() == 0
    }
    fn insert_nodes_for_node(
        &self,
        graph: &mut SimpleSeqGraph,
        node: NodeIndex,
        revcomp: bool,
    ) -> (NodeIndex, NodeIndex) {
        let weight = self.node_weight(node).unwrap();
        let seq = &weight.seq;
        let n = seq.len();
        let copy_num = weight.copy_num;
        let is_start_point_node = self.is_start_point_node(node);

        if revcomp {
            self.insert_nodes_for_seq(
                graph,
                &seq.to_revcomp(),
                |index, base| {
                    let pos = GenomeGraphPos::new(node, n - index - 1);
                    let is_start_point = is_start_point_node && index == 0;
                    SimpleSeqNode::new(copy_num, base, is_start_point, true, pos)
                },
                |_index| SimpleSeqEdge::new(Some(copy_num)),
            )
        } else {
            self.insert_nodes_for_seq(
                graph,
                seq,
                |index, base| {
                    let pos = GenomeGraphPos::new(node, index);
                    let is_start_point = is_start_point_node && index == 0;
                    SimpleSeqNode::new(copy_num, base, is_start_point, false, pos)
                },
                |_index| SimpleSeqEdge::new(Some(copy_num)),
            )
        }
    }
    ///
    /// convert a seq into nodes/edges in seq graph.
    ///
    /// seq ATCGATCG
    ///
    /// return tuple of NodeIndex of the first and last node.
    ///
    fn insert_nodes_for_seq<N, E>(
        &self,
        graph: &mut SimpleSeqGraph,
        seq: &Sequence,
        node_fn: N,
        edge_fn: E,
    ) -> (NodeIndex, NodeIndex)
    where
        N: Fn(usize, u8) -> SimpleSeqNode,
        E: Fn(usize) -> SimpleSeqEdge,
    {
        // add a new node for each bases
        let nodes: Vec<NodeIndex> = seq
            .iter()
            .enumerate()
            .map(|(index, &base)| graph.add_node(node_fn(index, base)))
            .collect();

        // add edges between adjacent two bases
        nodes
            .iter()
            .tuple_windows()
            .enumerate()
            .for_each(|(index, (&v0, &v1))| {
                graph.add_edge(v0, v1, edge_fn(index));
            });

        (*nodes.first().unwrap(), *nodes.last().unwrap())
    }
    ///
    /// Convert `GenomeGraph` into `SimpleSeqGraph`
    /// split a node containing bases into multiple nodes corresponding each bases
    ///
    pub fn to_seq_graph(&self) -> SimpleSeqGraph {
        self.to_seq_graph_with_revcomp(false)
    }
    ///
    /// Convert `GenomeGraph` into `SimpleSeqGraph`
    /// split a node containing bases into multiple nodes corresponding each bases
    ///
    pub fn to_seq_graph_with_revcomp(&self, with_revcomp: bool) -> SimpleSeqGraph {
        let mut graph = DiGraph::new();

        // mapping
        // from a node in genome graph (representing a sequence)
        // to head/tail of nodes in seq graph (each node representing a single base)
        //
        // mf: mapping of forward nodes
        // mb: mapping of backward (= revcomp) nodes
        //
        // Example
        // a genome graph node `ATCGG` will be converted into
        //
        // (forward)   A->T->C->G->G
        //             Head        Tail
        // (backward)  T<-A<-G<-C<-C
        //             Tail        Head
        let mut mf: HashMap<NodeIndex, (NodeIndex, NodeIndex)> = HashMap::new();
        let mut mb: HashMap<NodeIndex, (NodeIndex, NodeIndex)> = HashMap::new();

        // for each node
        for node in self.node_indices() {
            // forward
            let (head, tail) = self.insert_nodes_for_node(&mut graph, node, false);
            mf.insert(node, (head, tail));

            // backward
            if with_revcomp {
                let (head, tail) = self.insert_nodes_for_node(&mut graph, node, true);
                mb.insert(node, (head, tail));
            }
        }

        // for each edge
        for edge in self.edge_indices() {
            let (source, target) = self.edge_endpoints(edge).unwrap();
            let weight = self.edge_weight(edge).unwrap();

            // add edge for forward
            // (tail of source) -> (head of target)
            let (_, tail_of_source) = mf.get(&source).unwrap();
            let (head_of_target, _) = mf.get(&target).unwrap();
            graph.add_edge(
                *tail_of_source,
                *head_of_target,
                SimpleSeqEdge::new(weight.copy_num),
            );

            // add edge for backward
            // (tail of target) -> (head of source)
            if with_revcomp {
                let (_, tail_of_target) = mb.get(&target).unwrap();
                let (head_of_source, _) = mb.get(&source).unwrap();
                graph.add_edge(
                    *tail_of_target,
                    *head_of_source,
                    SimpleSeqEdge::new(weight.copy_num),
                );
            }
        }

        graph
    }
    /// Sample reads from the genome graph.
    pub fn sample_reads(&self, prof: &ReadProfile) -> Reads {
        self.sample_positioned_reads(prof).to_reads()
    }
    ///
    /// Sample positioned reads
    ///
    pub fn sample_positioned_reads(&self, prof: &ReadProfile) -> PositionedReads {
        // convert to phmm
        let sg = self.to_seq_graph();
        let phmm = sg.to_phmm(prof.phmm_params.clone());

        // determine automatically the starting node list
        // as a vector of node index in seqgraph.
        // and purge into prof.sample_profile.
        let mut prof = prof.clone();
        if let StartPoints::AllStartPoints = prof.sample_profile.start_points {
            let start_points = get_start_points(&sg);
            prof.sample_profile.start_points = StartPoints::Custom(start_points);
        }

        // sample reads using profile
        //
        // store the originated genome graph position in seqgraph
        let historys = phmm.sample_by_profile(&prof.sample_profile);
        historys.to_positioned_reads(&sg)
    }
    pub fn show_coverage(&self, reads: &PositionedReads) {
        let mut coverages: Vec<Vec<usize>> = self
            .0
            .node_indices()
            .map(|v| {
                let weight = self.0.node_weight(v).unwrap();
                let len = weight.seq.len();
                vec![0; len]
            })
            .collect();
        for read in reads.iter() {
            for origin in read.origins() {
                println!("{:?}", origin);
                coverages[origin.node.index()][origin.pos] += 1;
            }
        }

        for i in 0..self.node_count() {
            for j in 0..coverages[i].len() {
                println!("{} {} {}", i, j, coverages[i][j]);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ni, sequence_to_string};
    use crate::graph::seq_graph::find_node_from_source_pos;
    use crate::random_seq::generate;

    fn has_edge_between_two_pos(
        g: &DiGraph<SimpleSeqNode, SimpleSeqEdge>,
        source_a: GenomeGraphPos,
        source_b: GenomeGraphPos,
        is_revcomp: bool,
    ) -> bool {
        let a = find_node_from_source_pos(g, source_a, is_revcomp).unwrap();
        let b = find_node_from_source_pos(g, source_b, is_revcomp).unwrap();
        g.contains_edge(a, b)
    }

    #[test]
    fn genome_graph_linear() {
        let mut g = DiGraph::new();
        g.add_node(GenomeNode::new(b"GATGCTAGTT", 1));
        let gg = GenomeGraph(g);

        // forward only
        let sg = gg.to_seq_graph();
        assert_eq!(sg.node_count(), 10);
        assert_eq!(sg.edge_count(), 9);
        println!("{}", Dot::with_config(&sg, &[]));
        assert_eq!(get_start_points(&sg), vec![ni(0)]);

        // both strand
        let sg = gg.to_seq_graph_with_revcomp(true);
        println!("{}", Dot::with_config(&sg, &[]));
        let start_points = get_start_points(&sg);
        assert_eq!(start_points, vec![ni(0), ni(10)]);
        let start_point_sources: Vec<_> = start_points
            .iter()
            .map(|&v| sg.node_weight(v).unwrap().source())
            .collect();
        assert_eq!(
            start_point_sources,
            vec![GenomeGraphPos::new(ni(0), 0), GenomeGraphPos::new(ni(0), 9)]
        );
        assert_eq!(sg.node_count(), 20);
        assert_eq!(sg.edge_count(), 18);
    }

    #[test]
    fn genome_graph_branch() {
        let mut g = DiGraph::new();
        let v1 = g.add_node(GenomeNode::new(b"ATCGG", 1));
        let v2 = g.add_node(GenomeNode::new(b"AAC", 1));
        let v3 = g.add_node(GenomeNode::new(b"TTCG", 2));
        g.add_edge(v1, v3, GenomeEdge::new(None));
        g.add_edge(v2, v3, GenomeEdge::new(None));
        let gg = GenomeGraph(g);
        let sg = gg.to_seq_graph();
        assert_eq!(sg.node_count(), 12);
        assert_eq!(sg.edge_count(), 11);
        println!("{}", Dot::with_config(&sg, &[]));
        println!("{:?}", get_start_points(&sg));
        assert_eq!(get_start_points(&sg), vec![ni(0), ni(5)]);

        // both strand
        let sg = gg.to_seq_graph_with_revcomp(true);
        println!("{}", Dot::with_config(&sg, &[]));
        assert!(has_edge_between_two_pos(
            &sg,
            GenomeGraphPos::new(ni(0), 4),
            GenomeGraphPos::new(ni(2), 0),
            false
        ));
        assert!(has_edge_between_two_pos(
            &sg,
            GenomeGraphPos::new(ni(1), 2),
            GenomeGraphPos::new(ni(2), 0),
            false
        ));
        assert!(has_edge_between_two_pos(
            &sg,
            GenomeGraphPos::new(ni(2), 0),
            GenomeGraphPos::new(ni(0), 4),
            true
        ));
        assert!(has_edge_between_two_pos(
            &sg,
            GenomeGraphPos::new(ni(2), 0),
            GenomeGraphPos::new(ni(1), 2),
            true
        ));
    }

    #[test]
    fn genome_graph_circular() {
        let mut g = DiGraph::new();
        let v1 = g.add_node(GenomeNode::new(b"ATCGG", 1));
        g.add_edge(v1, v1, GenomeEdge::new(None));
        let gg = GenomeGraph(g);
        let sg = gg.to_seq_graph();
        assert_eq!(sg.node_count(), 5);
        assert_eq!(sg.edge_count(), 5);
        println!("{}", Dot::with_config(&sg, &[]));
        println!("{:?}", get_start_points(&sg));
        assert_eq!(get_start_points(&sg), vec![]);
    }

    #[test]
    fn genome_graph_from_seqs() {
        let seqs = vec![b"ATCGATTCGAT".to_vec(), b"CTCTTCTTCTCT".to_vec()];
        let gg = GenomeGraph::from_seqs(&seqs);
        assert_eq!(gg.node_count(), 2);
        assert_eq!(gg.edge_count(), 0);
    }

    #[test]
    fn genome_graph_sampling() {
        let seqs = vec![b"ATCGATTCGAT".to_vec(), b"CTCTTCTTCTCT".to_vec()];
        let graph = GenomeGraph::from_seqs(&seqs);
        let reads = graph.sample_reads(&ReadProfile {
            sample_profile: SampleProfile {
                read_amount: ReadAmount::Count(10),
                seed: 0,
                length: 1000,
                // start_points: StartPoints::Random,
                start_points: StartPoints::AllStartPoints,
                endable: false,
            },
            phmm_params: PHMMParams::default(),
        });
        for read in reads.iter() {
            println!("{}", sequence_to_string(read));
        }
    }
    #[test]
    fn genome_graph_sampling_random_genome() {
        let genome = vec![generate(500, 0)];
        let g = GenomeGraph::from_seqs(&genome);
        let profile = ReadProfile {
            sample_profile: SampleProfile {
                read_amount: ReadAmount::TotalBases(100),
                seed: 0,
                length: 100,
                start_points: StartPoints::Random,
                endable: false,
            },
            phmm_params: PHMMParams::default(),
        };
        let reads = g.sample_reads(&profile);
        assert_eq!(reads.len(), 2);
        for (i, read) in reads.iter().enumerate() {
            println!("{}", sequence_to_string(read));
            if i == 0 {
                assert_eq!(sequence_to_string(read), "TGAATCCTAGATCCCGTTGTCGGGGCTCGGCGTTTGCTTTCTTAGATTCCGATAAGTAGATGGTTTCCTGGGTGAGGGCACTATTAAAGCGGCGATTTG");
            } else if i == 1 {
                assert_eq!(sequence_to_string(read), "AGCGATTAAACACCCTATAAAAATGGCCATCCGCTGAGCTTGCATCACAGTTGGTCTTACACATGCCTGCTTCATCAAAGTCCCACTGCGCCATCA");
            }
        }
    }
}
