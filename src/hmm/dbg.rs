use super::base::{Node, PHMM};
use crate::dbg;
use crate::dbg::{DbgHash, DBG};
use crate::kmer::kmer::{tailing_kmers, Kmer};
use crate::prob::Prob;
use arrayvec::ArrayVec;

struct DbgPHMM {
    kmers: Vec<Kmer>,
    copy_nums: Vec<u32>,
    childs: Vec<ArrayVec<usize, 4>>,
    parents: Vec<ArrayVec<usize, 4>>,
    /*
    total_copy_num: u32,
    trans_probs: Vec<ArrayVec<Prob, 4>>,
    emissions: Vec<u8>,
    */
}

impl DbgPHMM {
    fn new(kmers: Vec<Kmer>, copy_nums: Vec<u32>) -> Option<DbgPHMM> {
        // construct DbgHash and
        // 1. check copy_num consistency
        // 2. add tailing kmers
        // 3. linearize kmers
        // let d = DbgHash::from(kmers, copy_nums);
        if !Self::is_copy_nums_valid(&kmers, &copy_nums) {
            return None;
        }
        let (childs, parents) = Self::create_childs_and_parents(&kmers);
        Self::add_tailing_kmers(&kmers, &parents);
        Some(DbgPHMM {
            kmers,
            copy_nums,
            childs,
            parents,
        })
    }
    fn create_childs_and_parents(
        kmers: &[Kmer],
    ) -> (Vec<ArrayVec<usize, 4>>, Vec<ArrayVec<usize, 4>>) {
        let mut childs = Vec::new();
        let mut parents = Vec::new();
        for kmer in kmers.iter() {
            let childs_of_kmer: ArrayVec<usize, 4> = kmers
                .iter()
                .enumerate()
                .filter(|(i, x)| x.prefix() == kmer.suffix())
                .map(|(i, x)| i)
                .collect();
            let parents_of_kmer: ArrayVec<usize, 4> = kmers
                .iter()
                .enumerate()
                .filter(|(i, x)| x.suffix() == kmer.prefix())
                .map(|(i, x)| i)
                .collect();
            println!(
                "kmer {} childs:{:?} parents:{:?}",
                kmer, childs_of_kmer, parents_of_kmer
            );
            childs.push(childs_of_kmer);
            parents.push(parents_of_kmer);
        }
        (childs, parents)
    }
    fn add_tailing_kmers(kmers: &[Kmer], parents: &[ArrayVec<usize, 4>]) {
        for (kmer, parents_of_kmer) in kmers.iter().zip(parents.iter()) {
            if parents_of_kmer.len() == 0 {
                // found no parent kmer
                let new_kmers = tailing_kmers(kmer);
                println!("add {:?}", new_kmers);
            }
        }
    }
    fn is_copy_nums_valid(kmers: &[Kmer], copy_nums: &[u32]) -> bool {
        // 1. check for copy_num validity
        for kmer in kmers.iter() {
            // in_copy
            let suffix = kmer.suffix();
            let in_copy: u32 = kmers
                .iter()
                .enumerate()
                .map(|(i, x)| {
                    if x.suffix() == suffix {
                        copy_nums[i]
                    } else {
                        0
                    }
                })
                .sum();
            // out_copy
            let out_copy: u32 = kmers
                .iter()
                .enumerate()
                .map(|(i, x)| {
                    if x.prefix() == suffix {
                        copy_nums[i]
                    } else {
                        0
                    }
                })
                .sum();
            if in_copy > 0 && out_copy > 0 && in_copy != out_copy {
                return false;
            }
        }
        true
    }
}
// impl PHMM for DbgPHMM {}

fn find_edges() {}

pub fn test() {
    println!("hello, i am dbg");
    /*
    let d = DbgPHMM::new(kmers, copy_nums).unwrap();

    let edge = Kmer::from(b"TCGA");
    let tailings = tailing_kmers(&edge);
    for t in tailings.iter() {
        println!("{}", t);
    }
    */
}
