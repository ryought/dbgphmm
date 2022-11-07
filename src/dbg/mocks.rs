use super::dbg::{Dbg, DbgEdge, DbgNode};
use super::hashdbg_v2::HashDbg;
use super::impls::{SimpleDbg, SimpleDbgEdge, SimpleDbgNode};
use crate::common::Sequence;
use crate::kmer::veckmer::{kmer, VecKmer};
use crate::random_seq::generate;
use itertools::Itertools;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

pub fn mock_base() -> SimpleDbg<VecKmer> {
    let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"ATCGGCT");
    SimpleDbg::from_hashdbg(&hd)
}

pub fn mock_simple() -> SimpleDbg<VecKmer> {
    let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"AAAGCTTGATT");
    SimpleDbg::from_hashdbg(&hd)
}

pub fn mock_two_seqs() -> SimpleDbg<VecKmer> {
    let mut hd: HashDbg<VecKmer> = HashDbg::new(4);
    hd.add_seq(b"AAAGCTTGATT");
    hd.add_seq(b"CGTATC");
    SimpleDbg::from_hashdbg(&hd)
}

pub fn mock_rep() -> SimpleDbg<VecKmer> {
    let mut hd: HashDbg<VecKmer> = HashDbg::new(4);
    hd.add_seq(b"AAAAAAAAAAAAA");
    hd.add_seq(b"CCCCCCCCCCCCCC");
    SimpleDbg::from_hashdbg(&hd)
}

///
/// AACTAGCTT x1
/// CCGTAGGGC x1
///
/// `TAG` is intersecting k-1mer.
///
pub fn mock_intersection() -> SimpleDbg<VecKmer> {
    let mut hd: HashDbg<VecKmer> = HashDbg::new(4);
    hd.add_seq(b"AACTAGCTT");
    hd.add_seq(b"CCGTAGGGC");
    SimpleDbg::from_hashdbg(&hd)
}

///
/// k=4 dbg
///
/// ATAGCT x1
/// TAAGCC x1
///
/// `TAGC` is intersecting.
///
pub fn mock_intersection_small() -> SimpleDbg<VecKmer> {
    let mut hd: HashDbg<VecKmer> = HashDbg::new(4);
    hd.add_seq(b"ATAGCT");
    hd.add_seq(b"TAAGCC");
    SimpleDbg::from_hashdbg(&hd)
}

///
/// crate a mock `SimpleDbg` from single random sequence of given length
///
pub fn mock_random_with_seq(k: usize, length: usize) -> (SimpleDbg<VecKmer>, Sequence) {
    let seq = generate(length, 1);
    let mut hd: HashDbg<VecKmer> = HashDbg::new(k);
    hd.add_seq(&seq);
    (SimpleDbg::from_hashdbg(&hd), seq)
}

pub fn mock_random(k: usize, length: usize) -> SimpleDbg<VecKmer> {
    let (dbg, _) = mock_random_with_seq(k, length);
    dbg
}

///
/// simple dbg from a manually constructed DiGraph
///
/// ATCGAGCATG
///
pub fn mock_manual() -> SimpleDbg<VecKmer> {
    let mut graph = DiGraph::new();

    // add nodes
    let nodes = vec![
        graph.add_node(SimpleDbgNode::new(kmer(b"nnnA"), 1)),
        graph.add_node(SimpleDbgNode::new(kmer(b"nnAT"), 1)),
        graph.add_node(SimpleDbgNode::new(kmer(b"nATC"), 1)),
        graph.add_node(SimpleDbgNode::new(kmer(b"ATCG"), 1)),
        graph.add_node(SimpleDbgNode::new(kmer(b"TCGA"), 1)),
        graph.add_node(SimpleDbgNode::new(kmer(b"CGAG"), 1)),
        graph.add_node(SimpleDbgNode::new(kmer(b"GAGC"), 1)),
        graph.add_node(SimpleDbgNode::new(kmer(b"AGCA"), 1)),
        graph.add_node(SimpleDbgNode::new(kmer(b"GCAT"), 1)),
        graph.add_node(SimpleDbgNode::new(kmer(b"CATG"), 1)),
        graph.add_node(SimpleDbgNode::new(kmer(b"ATGn"), 1)),
        graph.add_node(SimpleDbgNode::new(kmer(b"TGnn"), 1)),
        graph.add_node(SimpleDbgNode::new(kmer(b"Gnnn"), 1)),
    ];
    let n = nodes.len();

    // add edges
    let edges: Vec<_> = (0..n)
        .map(|i| graph.add_edge(nodes[i % n], nodes[(i + 1) % n], SimpleDbgEdge::new(None)))
        .collect();

    SimpleDbg::from_digraph(4, graph)
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dbg_mock_random() {
        let dbg = mock_random(8, 100);
        println!("{}", dbg);
    }
}
