//! copy_nums structs
//! copy_nums should be stored with reference to the cdbg
use crate::compressed_dbg::CompressedDBG;
use crate::graph::{Edgei, Node};
use crate::kmer::kmer::linear_seq_to_kmers;
use itertools::Itertools;
use log::warn;
use std::ops::{Index, IndexMut};

/// CopyNums stores
/// simple wrapper of copy_nums vector Vec<u32> etc
/// it allows index access by &Node
///
/// ```
/// use dbgphmm::copy_nums::CopyNums;
/// use dbgphmm::graph::Node;
/// // immutable
/// let c = CopyNums(vec![55, 22, 33]);
/// assert_eq!(c[&Node(1)], 22);
/// assert_eq!(c[&Node(2)], 33);
///
/// // mutable
/// let mut c2 = CopyNums(vec![55, 22, 33]);
/// c2[&Node(1)] = 2222;
/// assert_eq!(c2[&Node(1)], 2222);
/// ```
pub struct CopyNums(pub Vec<u32>);

impl Index<&Node> for CopyNums {
    type Output = u32;
    fn index(&self, v: &Node) -> &Self::Output {
        &self.0[v.0]
    }
}

impl IndexMut<&Node> for CopyNums {
    fn index_mut(&mut self, v: &Node) -> &mut Self::Output {
        &mut self.0[v.0]
    }
}

/// CopyNums constructors
impl CopyNums {
    pub fn zero(cdbg: &CompressedDBG) -> CopyNums {
        CopyNums(vec![0; cdbg.n_kmers()])
    }
    /// constructor from seqs
    pub fn from_seqs(cdbg: &CompressedDBG, seqs: &[Vec<u8>]) -> Option<CopyNums> {
        let mut c = CopyNums::zero(cdbg);
        for seq in seqs.iter() {
            for kmer in linear_seq_to_kmers(seq, cdbg.k()) {
                match cdbg.id(&kmer) {
                    Some(v) => c[&v] += 1,
                    None => {
                        warn!("true kmer {} not found in cdbg", kmer);
                        return None;
                    }
                }
            }
        }
        Some(c)
    }
}

/// EdgeCopyNums stores
/// simple wrapper of edge_copy_nums vector Vec<Vec<u32>> etc
/// it has same shape of cdbg.childs, so it can access by
/// &Edgei { source: Node, child_index: usize }
///
/// ```
/// use dbgphmm::copy_nums::EdgeCopyNums;
/// use dbgphmm::graph::{Node, Edgei};
/// let e = EdgeCopyNums(vec![vec![11, 33], vec![55], vec![99, 12]]);
/// assert_eq!(e[&Edgei::new(Node(0), 0)], 11);
/// assert_eq!(e[&Edgei::new(Node(0), 1)], 33);
/// assert_eq!(e[&Edgei::new(Node(2), 1)], 12);
///
/// let mut e2 = EdgeCopyNums(vec![vec![11, 33], vec![55], vec![99, 12]]);
/// e2[&Edgei::new(Node(0), 0)] = 11111;
/// assert_eq!(e2[&Edgei::new(Node(0), 0)], 11111);
/// ```
pub struct EdgeCopyNums(pub Vec<Vec<u32>>);

impl Index<&Edgei> for EdgeCopyNums {
    type Output = u32;
    fn index(&self, e: &Edgei) -> &Self::Output {
        &self.0[e.source.0][e.child_index]
    }
}

impl IndexMut<&Edgei> for EdgeCopyNums {
    fn index_mut(&mut self, e: &Edgei) -> &mut Self::Output {
        &mut self.0[e.source.0][e.child_index]
    }
}

/// EdgeCopyNums constructors
impl EdgeCopyNums {
    pub fn zero(cdbg: &CompressedDBG) -> EdgeCopyNums {
        let edge_copy_nums: Vec<Vec<u32>> = cdbg
            .iter_nodes()
            .map(|v| vec![0; cdbg.childs(&v).len()])
            .collect();
        EdgeCopyNums(edge_copy_nums)
    }
    /// constructor from seqs
    pub fn from_seqs(cdbg: &CompressedDBG, seqs: &[Vec<u8>]) {
        // let mut copy_nums: Vec<u32> = vec![0; cdbg.n_kmers()];
        for seq in seqs.iter() {
            for (k0, k1) in linear_seq_to_kmers(seq, cdbg.k()).tuple_windows() {
                match (cdbg.id(&k0), cdbg.id(&k1)) {
                    (Some(v0), Some(v1)) => {
                        println!("k0={}({:?}) k1={}({:?})", k0, v0, k1, v1);
                    }
                    _ => {
                        warn!("true kmer {} or {} not found in cdbg", k0, k1);
                        // return None;
                    }
                }
            }
        }
        // Some(CopyNums(copy_nums))
    }
}

// node copy-num labeled
pub struct Ncdbg<'a> {
    pub cdbg: &'a CompressedDBG,
    pub copy_nums: CopyNums,
}

impl<'a> Ncdbg<'a> {
    pub fn new(cdbg: &CompressedDBG, copy_nums: CopyNums) -> Ncdbg {
        assert!(Ncdbg::is_consistent(cdbg, &copy_nums));
        Ncdbg { cdbg, copy_nums }
    }
    fn is_consistent(cdbg: &CompressedDBG, copy_nums: &CopyNums) -> bool {
        cdbg.is_consistent_copy_num(&copy_nums.0)
    }
    pub fn copy_num(&self, v: &Node) -> u32 {
        if !self.cdbg.is_valid_node(v) {
            panic!("node {:?} is not in cdbg", v);
        }
        self.copy_nums[v]
    }
}

// edge copy-num labeled
pub struct Ecdbg<'a> {
    pub cdbg: &'a CompressedDBG,
    pub node_copy_nums: CopyNums,
    pub edge_copy_nums: EdgeCopyNums,
}

impl<'a> Ecdbg<'a> {
    pub fn new(
        cdbg: &CompressedDBG,
        node_copy_nums: CopyNums,
        edge_copy_nums: EdgeCopyNums,
    ) -> Ecdbg {
        Ecdbg {
            cdbg,
            node_copy_nums,
            edge_copy_nums,
        }
    }
    pub fn copy_num(&self, v: &Node) -> u32 {
        if !self.cdbg.is_valid_node(v) {
            panic!("node {:?} is not in cdbg", v);
        }
        self.node_copy_nums[v]
    }
    pub fn edge_copy_num(&self, v: &Node, w: &Node) -> u32 {
        if !self.cdbg.is_adjacent(v, w) {
            panic!("edge (v,w) = ({:?},{:?}) is not in cdbg", v, w);
        } else {
            let i = self.cdbg.childs(v).iter().position(|u| u == w).unwrap();
            self.edge_copy_nums[&Edgei::new(*v, i)]
        }
    }
    fn is_node_consistent(&self) -> bool {
        self.cdbg.is_consistent_copy_num(&self.node_copy_nums.0)
    }
    /// check that edge_{child,parent}_consistency
    /// for all nodes
    fn is_edge_consistent(&self) -> bool {
        for v in self.cdbg.iter_nodes() {
            if !self.is_edge_child_consistent(&v) || !self.is_edge_parent_consistent(&v) {
                return false;
            }
        }
        true
    }
    fn is_edge_child_consistent(&self, v: &Node) -> bool {
        let sum_copy_number_of_childs: u32 = self
            .cdbg
            .childs(v)
            .iter()
            .map(|w| self.edge_copy_num(v, w))
            .sum();
        let my_copy_number: u32 = self.copy_num(v);
        my_copy_number == sum_copy_number_of_childs
    }
    fn is_edge_parent_consistent(&self, v: &Node) -> bool {
        let sum_copy_number_of_parents: u32 = self
            .cdbg
            .parents(v)
            .iter()
            .map(|w| self.edge_copy_num(w, v))
            .sum();
        let my_copy_number: u32 = self.copy_num(v);
        my_copy_number == sum_copy_number_of_parents
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mocks::test_cdbg_01;
    #[test]
    fn ncdbg_0() {
        let (cdbg, cn, _) = test_cdbg_01();
        let ncdbg = Ncdbg::new(&cdbg, cn);
        println!("{}", ncdbg.copy_num(&Node(0)));
    }
    #[test]
    fn ecdbg_0() {
        let (cdbg, cn, ecn) = test_cdbg_01();
        let ecdbg = Ecdbg::new(&cdbg, cn, ecn);
        assert!(ecdbg.is_node_consistent());
        assert!(ecdbg.is_edge_consistent());
    }
    #[test]
    fn edgecopynums_fromseqs() {
        let seqs: Vec<Vec<u8>> = vec![b"ATTCGATCGATTT".to_vec()];
        let (cdbg, cn, ecn) = test_cdbg_01();
        EdgeCopyNums::from_seqs(&cdbg, &seqs);
    }
}
