//! copy_nums structs
//! copy_nums should be stored with reference to the cdbg
use crate::compressed_dbg::CompressedDBG;
use crate::graph::{Edgep, Node};
use std::ops::Index;

pub struct CopyNums(pub Vec<u32>);
pub struct EdgeCopyNums(pub Vec<Vec<u32>>);
impl Index<&Node> for CopyNums {
    type Output = u32;
    fn index(&self, v: &Node) -> &Self::Output {
        &self.0[v.0]
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
        // TODO assert v is valid node
        self.copy_nums[v]
    }
}

// edge copy-num labeled
pub struct Ecdbg<'a> {
    pub cdbg: &'a CompressedDBG,
    pub node_copy_nums: Vec<u32>,
    pub edge_copy_nums: Vec<Vec<u32>>,
}

impl<'a> Ecdbg<'a> {
    pub fn new(
        cdbg: &CompressedDBG,
        node_copy_nums: Vec<u32>,
        edge_copy_nums: Vec<Vec<u32>>,
    ) -> Ecdbg {
        assert!(Ecdbg::is_node_consistent(cdbg, &node_copy_nums));
        assert!(Ecdbg::is_edge_consistent(
            cdbg,
            &node_copy_nums,
            &edge_copy_nums
        ));
        Ecdbg {
            cdbg,
            node_copy_nums,
            edge_copy_nums,
        }
    }
    fn is_node_consistent(cdbg: &CompressedDBG, node_copy_nums: &[u32]) -> bool {
        cdbg.is_consistent_copy_num(node_copy_nums)
    }
    fn is_edge_consistent(_: &CompressedDBG, _: &[u32], _: &[Vec<u32>]) -> bool {
        // TODO
        true
    }
    fn copy_num(&self, v: &Node) -> u32 {
        self.node_copy_nums[v.0]
    }
    /*
    fn edge_copy_num(&self, v: &Node, w: &Node) -> u32 {
        if !self.cdbg.is_adjacent(v, w) {
            panic!("edge (v,w) = ({:?},{:?}) is not in cdbg", v, w);
        } else {
            let i = self.cdbg.childs(v).iter().position(|u| u == w).unwrap();
            self.edge_copy_nums[]
        }
        0
    }
    */
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mocks::test_cdbg_01;
    #[test]
    fn ncdbg_0() {
        let (cdbg, copy_nums) = test_cdbg_01();
        println!("n={}, c={}", cdbg.n_kmers(), cdbg.n_cycles());
        let ncdbg = Ncdbg::new(&cdbg, CopyNums(copy_nums));
        println!("{}", ncdbg.copy_num(&Node(0)));
    }
}
