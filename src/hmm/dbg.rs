use super::base::{Node, PHMM};
use super::params::PHMMParams;
use crate::dbg;
use crate::dbg::{DbgHash, DBG};
use crate::kmer::kmer::{tailing_kmers, Kmer};
use crate::prob::Prob;
use arrayvec::ArrayVec;

pub struct DbgPHMM {
    dbg: DbgHash,
    // from vectorize
    kmers: Vec<Kmer>,
    childs: Vec<Vec<usize>>,
    parents: Vec<Vec<usize>>,
    trans_probs: Vec<Vec<Prob>>,
    // add manually
    copy_nums: Vec<u32>,
    emissions: Vec<u8>,
}

impl DbgPHMM {
    pub fn new(kmers: Vec<Kmer>, copy_nums: Vec<u32>) -> Option<DbgPHMM> {
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
            dbg: d,
            kmers,
            childs,
            parents,
            trans_probs,
            copy_nums,
            emissions,
        })
    }
}

impl PHMM for DbgPHMM {
    fn nodes(&self) -> Vec<Node> {
        (0..self.kmers.len()).map(|i| Node(i)).collect()
    }
    fn childs(&self, v: &Node) -> Vec<Node> {
        self.childs[v.0].iter().map(|i| Node(*i)).collect()
    }
    fn parents(&self, v: &Node) -> Vec<Node> {
        self.parents[v.0].iter().map(|i| Node(*i)).collect()
    }
    fn is_adjacent(&self, v: &Node, w: &Node) -> bool {
        self.childs[v.0].iter().any(|&i| i == w.0)
    }
    fn copy_num(&self, v: &Node) -> u32 {
        self.copy_nums[v.0]
    }
    fn emission(&self, v: &Node) -> u8 {
        self.emissions[v.0]
    }
    fn trans_prob(&self, v: &Node, w: &Node) -> Prob {
        match self.childs[v.0].iter().position(|&i| i == w.0) {
            Some(index) => self.trans_probs[v.0][index],
            None => Prob::from_prob(0.0),
        }
    }
}

pub fn test() {
    let kmers: Vec<Kmer> = vec![
        Kmer::from(b"GGAC"),
        Kmer::from(b"TGAC"),
        Kmer::from(b"GACT"),
        Kmer::from(b"GACC"),
        Kmer::from(b"ACCT"),
        Kmer::from(b"CCTG"),
    ];
    let copy_nums: Vec<u32> = vec![1, 2, 2, 1, 1, 1];
    let d = DbgPHMM::new(kmers, copy_nums).unwrap();
    // println!("{}", d.dbg.as_dot());
    // println!("{}", d.as_dot());

    /*
    let param = PHMMParams::new(
        Prob::from_prob(0.01),
        Prob::from_prob(0.01),
        Prob::from_prob(0.01),
        10,
    );
    let q = b"TGACCT";
    let p = d.forward_prob(&param, q);
    println!("prob = {}", p);
    */
}
