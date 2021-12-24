use crate::compressed_dbg::CompressedDBG;
use crate::copy_nums::{CopyNums, EdgeCopyNums};
use crate::dbg::{DbgHash, DBG};

/// k=8 mock cdbg (and its copy nums) representing `ATTCGATTCGATCGGA`
/// with no loops
pub fn test_cdbg_01() -> (CompressedDBG, CopyNums, EdgeCopyNums) {
    let seqs: Vec<Vec<u8>> = vec![b"ATTCGATCGATTT".to_vec()];
    let (cdbg, cn) = CompressedDBG::from_seqs(&seqs, 8);
    // TODO separate this creation for random some valid edge weights
    let ecn: Vec<Vec<u32>> = cdbg
        .iter_nodes()
        .map(|v| cdbg.childs(&v).iter().map(|_| 1).collect())
        .collect();
    (cdbg, CopyNums(cn), EdgeCopyNums(ecn))
}

/// k=8 mock cdbg (and its copy nums) representing `ATCGATTCGATCGATTCGATAGATCG`
/// with 2 cycles
pub fn test_cdbg_02() -> (CompressedDBG, Vec<u32>) {
    let seqs: Vec<Vec<u8>> = vec![b"ATCGATTCGATCGATTCGATAGATCG".to_vec()];
    CompressedDBG::from_seqs(&seqs, 8)
}
