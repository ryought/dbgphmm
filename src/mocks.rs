use crate::compressed_dbg::CompressedDBG;
use crate::dbg::{DbgHash, DBG};

/// k=8 mock cdbg (and its copy nums) representing `ATTCGATTCGATCGGA`
/// with no loops
pub fn test_cdbg_01() -> (CompressedDBG, Vec<u32>) {
    let seqs: Vec<Vec<u8>> = vec![b"ATTCGATCGATTT".to_vec()];
    CompressedDBG::from_seqs(&seqs, 8)
}

/// k=8 mock cdbg (and its copy nums) representing `ATCGATTCGATCGATTCGATAGATCG`
/// with 2 cycles
pub fn test_cdbg_02() -> (CompressedDBG, Vec<u32>) {
    let seqs: Vec<Vec<u8>> = vec![b"ATCGATTCGATCGATTCGATAGATCG".to_vec()];
    CompressedDBG::from_seqs(&seqs, 8)
}
