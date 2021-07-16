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

pub struct FCDbgPHMM<'a> {
    cdbg: &'a CompressedDBG,
    freqs: Vec<f64>,
    total_freq: f64,
    trans_probs: Vec<Vec<Prob>>,
}

impl<'a> FCDbgPHMM<'a> {
    /// construct FCDbgPHMM from frequencies
    pub fn new(cdbg: &CompressedDBG, freqs: Vec<f64>) -> FCDbgPHMM {
        FCDbgPHMM {
            cdbg,
            trans_probs: cdbg.freq_to_trans_prob(&freqs),
            total_freq: cdbg.total_emitable_freq(&freqs),
            freqs,
        }
    }
}

impl<'a> PHMM for FCDbgPHMM<'a> {
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
        self.freqs[v.0]
    }
    fn total_copy_num(&self) -> f64 {
        self.total_freq
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

impl<'a> PHMMSampler for FCDbgPHMM<'a> {}
