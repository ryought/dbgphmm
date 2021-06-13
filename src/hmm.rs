pub mod base;
pub mod linear;
pub mod testing;

/*
struct DbgPHMM {
    kmers: Vec<Kmer>,
    copy_nums: Vec<u32>,
    total_copy_num: u32,
    childs: Vec<ArrayVec>,
    parents: Vec<ArrayVec>,
    trans_probs: Vec<ArrayVec>,
    emissions: Vec<u8>,
}
impl DbgPHMM {
    pub fn from(kmers, copy_nums) -> DbgPHMM {
        // 1. calc flow
        // 2. store adjacency
    }
}
impl PHMM for DbgPHMM {}
*/
