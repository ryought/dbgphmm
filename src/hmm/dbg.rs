use super::base::{Node, PHMM};
use super::params::PHMMParams;
use super::sampler::PHMMSampler;
use crate::dbg;
use crate::dbg::{DbgHash, DBG};
use crate::kmer::kmer::Kmer;
use crate::prob::Prob;
use arrayvec::ArrayVec;
use log::{info, warn};

pub struct DbgPHMM {
    pub dbg: DbgHash,
    // from vectorize
    kmers: Vec<Kmer>,
    nodes: Vec<Node>,
    childs: Vec<ArrayVec<Node, 5>>,
    parents: Vec<ArrayVec<Node, 5>>,
    trans_probs: Vec<ArrayVec<Prob, 5>>,
    // add manually
    copy_nums: Vec<u32>,
    total_copy_num: u32,
    emissions: Vec<u8>,
}

impl DbgPHMM {
    pub fn from_dbg(dbg: DbgHash) -> DbgPHMM {
        // linearize kmers
        let (kmers, childs, parents, trans_probs) = dbg.vectorize();
        let copy_nums: Vec<u32> = kmers.iter().map(|kmer| dbg.find(kmer)).collect();
        let total_copy_num: u32 = copy_nums.iter().sum();
        let emissions: Vec<u8> = kmers.iter().map(|kmer| kmer.last()).collect();
        let nodes: Vec<Node> = (0..kmers.len()).map(|i| Node(i)).collect();

        DbgPHMM {
            dbg,
            kmers,
            nodes,
            childs,
            parents,
            trans_probs,
            copy_nums,
            total_copy_num,
            emissions,
        }
    }
    pub fn from_seqs(seqs: Vec<Vec<u8>>, k: usize) -> DbgPHMM {
        // construct dbg
        let d = DbgHash::from_seqs(seqs, k);
        DbgPHMM::from_dbg(d)
    }
    pub fn new(kmers: Vec<Kmer>, copy_nums: Vec<u32>) -> Option<DbgPHMM> {
        // construct DbgHash and
        let mut d = DbgHash::from(kmers, copy_nums);

        // 1. check copy_num consistency
        if !d.is_copy_number_consistent() {
            // return None;
        }

        // 2. add tailing kmers
        d.augment_edge_kmers();

        // 3. linearize kmers
        Some(DbgPHMM::from_dbg(d))
    }
}

impl PHMM for DbgPHMM {
    fn n_nodes(&self) -> usize {
        self.nodes.len()
    }
    fn childs(&self, v: &Node) -> &[Node] {
        &self.childs[v.0]
    }
    fn parents(&self, v: &Node) -> &[Node] {
        &self.parents[v.0]
    }
    fn is_adjacent(&self, v: &Node, w: &Node) -> bool {
        self.childs[v.0].iter().any(|&i| i == *w)
    }
    fn copy_num(&self, v: &Node) -> u32 {
        self.copy_nums[v.0]
    }
    fn total_copy_num(&self) -> u32 {
        self.total_copy_num
    }
    fn emission(&self, v: &Node) -> u8 {
        self.emissions[v.0]
    }
    fn trans_prob(&self, v: &Node, w: &Node) -> Prob {
        match self.childs[v.0].iter().position(|&i| i == *w) {
            Some(index) => self.trans_probs[v.0][index],
            None => Prob::from_prob(0.0),
        }
    }
}

impl PHMMSampler for DbgPHMM {}

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
