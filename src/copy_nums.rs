use crate::compressed_dbg::CompressedDBG;

// node copy-num labeled
pub struct Ncdbg<'a> {
    pub cdbg: &'a CompressedDBG,
    pub copy_nums: Vec<u32>,
}

impl<'a> Ncdbg<'a> {
    pub fn new(cdbg: &CompressedDBG, copy_nums: Vec<u32>) -> Ncdbg {
        assert!(Ncdbg::is_consistent(cdbg, &copy_nums));
        Ncdbg { cdbg, copy_nums }
    }
    fn is_consistent(cdbg: &CompressedDBG, copy_nums: &[u32]) -> bool {
        cdbg.is_consistent_copy_num(copy_nums)
    }
}

// edge copy-num labeled
pub struct Ecdbg<'a> {
    pub cdbg: &'a CompressedDBG,
    pub node_copy_nums: Vec<u32>,
    pub edge_copy_nums: Vec<u32>,
}

impl<'a> Ecdbg<'a> {
    pub fn new(cdbg: &CompressedDBG, node_copy_nums: Vec<u32>, edge_copy_nums: Vec<u32>) -> Ecdbg {
        assert!(Ecdbg::is_node_consistent(cdbg, &node_copy_nums));
        assert!(Ecdbg::is_edge_consistent());
        Ecdbg {
            cdbg,
            node_copy_nums,
            edge_copy_nums,
        }
    }
    fn is_node_consistent(cdbg: &CompressedDBG, node_copy_nums: &[u32]) -> bool {
        cdbg.is_consistent_copy_num(node_copy_nums)
    }
    fn is_edge_consistent() -> bool {
        // TODO
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mocks::test_cdbg_01;
    #[test]
    fn ncdbg_0() {
        let (cdbg, copy_nums) = test_cdbg_01();
        println!("n={}, c={}", cdbg.n_kmers(), cdbg.n_cycles());
        let ncdbg = Ncdbg::new(&cdbg, copy_nums);
    }
}
