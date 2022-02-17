use super::{EDbg, EDbgEdge, EDbgNode};
use crate::common::{CopyNum, Sequence};
use crate::kmer::kmer::{Kmer, KmerLike};
use petgraph::graph::DiGraph;

/// Basic implementations of EDbg
pub type SimpleEDbg<K> = EDbg<SimpleEDbgNode, SimpleEDbgEdge<K>>;

/// Basic implementations of EDbgNode
pub struct SimpleEDbgNode();

impl EDbgNode for SimpleEDbgNode {}

impl SimpleEDbgNode {
    pub fn new() -> Self {
        SimpleEDbgNode()
    }
}

/// Basic implementations of EDbgNode
pub struct SimpleEDbgEdge<K: KmerLike> {
    kmer: K,
    copy_num: CopyNum,
}

impl<K: KmerLike> EDbgEdge for SimpleEDbgEdge<K> {
    type Kmer = K;
    fn kmer(&self) -> &K {
        &self.kmer
    }
    fn copy_num(&self) -> CopyNum {
        self.copy_num
    }
}

impl<K: KmerLike> SimpleEDbgEdge<K> {
    pub fn new(kmer: K, copy_num: CopyNum) -> Self {
        SimpleEDbgEdge { kmer, copy_num }
    }
}

impl std::fmt::Display for SimpleEDbgNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "")
    }
}

impl<K: KmerLike> std::fmt::Display for SimpleEDbgEdge<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} (x{})", self.kmer, self.copy_num)
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::veckmer::VecKmer;

    #[test]
    fn edbg_simple() {
        let mut graph = DiGraph::new();
        let v1 = graph.add_node(SimpleEDbgNode::new());
        let v2 = graph.add_node(SimpleEDbgNode::new());
        let e1 = graph.add_edge(v1, v2, SimpleEDbgEdge::new(VecKmer::from_bases(b"ATCG"), 1));
        let e1 = graph.add_edge(v1, v2, SimpleEDbgEdge::new(VecKmer::from_bases(b"ATCG"), 1));
        let dbg: SimpleEDbg<VecKmer> = SimpleEDbg::new(4, graph);
    }
}
