use super::base::{Node, PHMM};
use super::sampler::PHMMSampler;
use crate::prob::Prob;

pub struct LinearPHMM {
    bases: Vec<u8>,
}
impl LinearPHMM {
    pub fn from(seq: &[u8]) -> LinearPHMM {
        LinearPHMM {
            bases: seq.to_vec(),
        }
    }
}
impl PHMM for LinearPHMM {
    fn nodes(&self) -> Vec<Node> {
        (0..self.bases.len()).map(|i| Node(i)).collect()
    }
    fn childs(&self, v: &Node) -> Vec<Node> {
        if v.0 != self.nodes().len() - 1 {
            vec![Node(v.0 + 1)]
        } else {
            vec![]
        }
    }
    fn parents(&self, v: &Node) -> Vec<Node> {
        if v.0 != 0 {
            vec![Node(v.0 - 1)]
        } else {
            vec![]
        }
    }
    fn is_adjacent(&self, v: &Node, w: &Node) -> bool {
        v.0 + 1 == w.0
    }
    fn copy_num(&self, v: &Node) -> u32 {
        if v.0 < self.bases.len() {
            1
        } else {
            0
        }
    }
    fn emission(&self, v: &Node) -> u8 {
        self.bases[v.0]
    }
    fn trans_prob(&self, v: &Node, w: &Node) -> Prob {
        if self.is_adjacent(v, w) {
            Prob::from_prob(1.0)
        } else {
            Prob::from_prob(0.0)
        }
    }
}

impl PHMMSampler for LinearPHMM {}
