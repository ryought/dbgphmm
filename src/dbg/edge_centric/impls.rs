use super::{EDbg, EDbgEdge, EDbgNode};
use crate::common::{CopyNum, Sequence};
use crate::kmer::kmer::{Kmer, KmerLike};

/// Basic implementations of EDbg
pub type SimpleEDbg<K> = EDbg<SimpleEDbgNode<K>, SimpleEDbgEdge<K>>;

/// Basic implementations of EDbgNode
pub struct SimpleEDbgNode<K: KmerLike> {
    km1mer: K,
    copy_num: CopyNum,
}

impl<K: KmerLike> EDbgNode for SimpleEDbgNode<K> {
    type Kmer = K;
    fn km1mer(&self) -> &K {
        &self.km1mer
    }
    fn copy_num(&self) -> CopyNum {
        self.copy_num
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

impl<K: KmerLike> std::fmt::Display for SimpleEDbgNode<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} (x{})", self.km1mer, self.copy_num)
    }
}

impl<K: KmerLike> std::fmt::Display for SimpleEDbgEdge<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} (x{})", self.kmer, self.copy_num)
    }
}
