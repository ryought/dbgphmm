use crate::compressed_dbg::CompressedDBG;

// node copy-num labeled
struct Ncdbg<'a> {
    cdbg: &'a CompressedDBG,
    copy_nums: Vec<u32>,
}

impl<'a> Ncdbg<'a> {
    pub fn new(cdbg: &CompressedDBG, copy_nums: Vec<u32>) -> Ncdbg {
        Ncdbg { cdbg, copy_nums }
    }
}

// edge copy-num labeled
struct Ecdbg<'a> {
    cdbg: &'a CompressedDBG,
    node_copy_nums: Vec<u32>,
    edge_copy_nums: Vec<u32>,
}

impl<'a> Ecdbg<'a> {
    pub fn new(cdbg: &CompressedDBG, node_copy_nums: Vec<u32>, edge_copy_nums: Vec<u32>) -> Ecdbg {
        // TODO assert edge copy num is consistent to node copy num
        Ecdbg {
            cdbg,
            node_copy_nums,
            edge_copy_nums,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mocks::test_cdbg_01;
    #[test]
    fn ncdbg_0() {
        let (cdbg, _) = test_cdbg_01();
        println!("n={}, c={}", cdbg.n_kmers(), cdbg.n_cycles())
    }
}
