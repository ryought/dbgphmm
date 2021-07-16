use super::base::{Node, PHMM};
use super::params::PHMMParams;
use super::sampler::PHMMSampler;
use crate::compressed_dbg::CompressedDBG;
use crate::dbg;
use crate::dbg::{DbgHash, DBG};
use crate::kmer::kmer::Kmer;
use crate::prob::Prob;
use arrayvec::ArrayVec;
use log::{info, warn};

pub struct CDbgPHMM<'a> {
    cdbg: &'a CompressedDBG,
    copy_nums: Vec<u32>,
    total_copy_num: u32,
    trans_probs: Vec<Vec<Prob>>,
}

impl<'a> CDbgPHMM<'a> {
    /// construct CDbgPHMM from cdbg and their copy_nums
    pub fn new(cdbg: &CompressedDBG, copy_nums: Vec<u32>) -> CDbgPHMM {
        let trans_probs = cdbg.copy_num_to_trans_prob(&copy_nums);
        // TODO this should omit the emission less nodes or not?
        let total_copy_num = cdbg.total_emitable_copy_num(&copy_nums);
        CDbgPHMM {
            cdbg,
            copy_nums,
            trans_probs,
            total_copy_num,
        }
    }
    /// check if CDbgPHMM copy numbers (flows) are consistent or not
    /// FIXME to be implemented
    pub fn is_consistent(&self) -> bool {
        true
    }
}

impl<'a> PHMM for CDbgPHMM<'a> {
    fn n_nodes(&self) -> usize {
        self.cdbg.n_kmers()
    }
    fn childs(&self, v: &Node) -> &[Node] {
        self.cdbg.childs(v)
    }
    fn parents(&self, v: &Node) -> &[Node] {
        self.cdbg.parents(v)
    }
    fn is_adjacent(&self, v: &Node, w: &Node) -> bool {
        self.cdbg.childs(v).iter().any(|&i| i == *w)
    }
    fn copy_num(&self, v: &Node) -> f64 {
        self.copy_nums[v.0] as f64
    }
    fn total_copy_num(&self) -> f64 {
        self.total_copy_num as f64
    }
    fn emission(&self, v: &Node) -> u8 {
        self.cdbg.emission(v)
    }
    fn trans_prob(&self, v: &Node, w: &Node) -> Prob {
        match self.cdbg.childs(v).iter().position(|&i| i == *w) {
            Some(index) => self.trans_probs[v.0][index],
            None => Prob::from_prob(0.0),
        }
    }
}

impl<'a> PHMMSampler for CDbgPHMM<'a> {}

pub fn test() {
    let mut d = DbgHash::new();
    // d.add_seq(b"ATCGATTCGATCGATTCGATAGATCG", 8);
    d.add_seq(b"ATCGATTCGATCG", 8);
    let cdbg = CompressedDBG::from(&d, 8);
    let copy_nums_true: Vec<u32> = cdbg
        .iter_nodes()
        .map(|v| {
            let kmer = cdbg.kmer(&v);
            d.find(kmer)
        })
        .collect();

    let param = PHMMParams::default();
    let read = b"GGATCGA";

    // 0. construct cdbgphmm
    let model = CDbgPHMM::new(&cdbg, copy_nums_true);
    // 1. check graphviz
    // println!("{}", model.as_dot());
    println!("{}", model.forward_prob(&param, read));

    // 2. check prob is same
    let old_model = super::dbg::DbgPHMM::from_dbg(d);
    // println!("{}", old_model.as_dot());
    println!("{}", old_model.forward_prob(&param, read));
}
