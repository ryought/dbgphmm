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
    fn copy_num(&self, v: &Node) -> u32 {
        self.copy_nums[v.0]
    }
    fn total_copy_num(&self) -> u32 {
        self.total_copy_num
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
