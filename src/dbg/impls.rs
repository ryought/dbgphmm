use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::common::{CopyNum, Sequence};
use crate::kmer::kmer::{Kmer, KmerLike};

/// Basic implementations of Dbg
pub type SimpleDbg<K> = Dbg<SimpleDbgNode<K>, SimpleDbgEdge>;

/// Basic implementations of DbgNode
pub struct SimpleDbgNode<K: KmerLike> {
    kmer: K,
    copy_num: CopyNum,
}

impl<K: KmerLike> DbgNode for SimpleDbgNode<K> {
    type Kmer = K;
    fn new(kmer: K, copy_num: CopyNum) -> Self {
        SimpleDbgNode { kmer, copy_num }
    }
    fn kmer(&self) -> &K {
        &self.kmer
    }
    fn copy_num(&self) -> CopyNum {
        self.copy_num
    }
    fn set_copy_num(&mut self, copy_num: CopyNum) {
        self.copy_num = copy_num;
    }
}

/// Basic implementations of DbgNode
pub struct SimpleDbgEdge {
    copy_num: Option<CopyNum>,
}

impl DbgEdge for SimpleDbgEdge {
    fn new(copy_num: Option<CopyNum>) -> Self {
        SimpleDbgEdge { copy_num }
    }
    fn copy_num(&self) -> Option<CopyNum> {
        self.copy_num
    }
    fn set_copy_num(&mut self, copy_num: Option<CopyNum>) {
        self.copy_num = copy_num;
    }
}

impl<K: KmerLike> std::fmt::Display for SimpleDbgNode<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} (x{})", self.kmer(), self.copy_num)
    }
}

impl std::fmt::Display for SimpleDbgEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.copy_num {
            Some(copy_num) => write!(f, "x{}", copy_num),
            None => write!(f, "?"),
        }
    }
}
