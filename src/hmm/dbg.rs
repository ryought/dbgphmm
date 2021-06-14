use super::base::{Node, PHMM};
use crate::dbg;
use crate::dbg::{DbgHash, DBG};
use crate::kmer::kmer::{tailing_kmers, Kmer};
use crate::prob::Prob;
use arrayvec::ArrayVec;

struct DbgPHMM {
    // from vectorize
    kmers: Vec<Kmer>,
    childs: Vec<Vec<usize>>,
    parents: Vec<Vec<usize>>,
    trans_probs: Vec<Vec<Prob>>,
    // add manually
    copy_nums: Vec<u32>,
    emissions: Vec<u8>,
    // total_copy_num: u32,
}

impl DbgPHMM {
    fn new(kmers: Vec<Kmer>, copy_nums: Vec<u32>) -> Option<DbgPHMM> {
        // construct DbgHash and
        let mut d = DbgHash::from(kmers, copy_nums);

        // 1. check copy_num consistency
        if !d.is_copy_number_consistent() {
            return None;
        }

        // 2. add tailing kmers
        d.augment_edge_kmers();

        // 3. linearize kmers
        let (kmers, childs, parents, trans_probs) = d.vectorize();
        let copy_nums: Vec<u32> = kmers.iter().map(|kmer| d.find(kmer)).collect();
        let emissions: Vec<u8> = kmers.iter().map(|kmer| kmer.last()).collect();

        Some(DbgPHMM {
            kmers,
            childs,
            parents,
            trans_probs,
            copy_nums,
            emissions,
        })
    }
}
// impl PHMM for DbgPHMM {}

pub fn test() {
    println!("hello, i am dbg");
}
