use super::base::{Node, PHMM};
use super::sampler::PHMMSampler;
use crate::prob::Prob;

pub struct LinearPHMM {
    bases: Vec<u8>,
    nodes: Vec<Node>,
}
impl LinearPHMM {
    pub fn from(seq: &[u8]) -> LinearPHMM {
        LinearPHMM {
            bases: seq.to_vec(),
            nodes: (0..seq.len()).map(|i| Node(i)).collect(),
        }
    }
}
impl PHMM for LinearPHMM {
    fn n_nodes(&self) -> usize {
        self.nodes.len()
    }
    fn childs(&self, v: &Node) -> &[Node] {
        // let mut childs = ArrayVec::<Node, 4>::new();
        if v.0 == self.nodes.len() - 1 {
            &[]
        } else {
            &self.nodes[v.0 + 1..v.0 + 2]
        }
    }
    fn parents(&self, v: &Node) -> &[Node] {
        // let mut childs = ArrayVec::<Node, 4>::new();
        if v.0 == 0 {
            &[]
        } else {
            &self.nodes[v.0 - 1..v.0]
        }
    }
    fn is_adjacent(&self, v: &Node, w: &Node) -> bool {
        v.0 + 1 == w.0
    }
    fn copy_num(&self, v: &Node) -> f64 {
        if v.0 < self.bases.len() {
            1.0
        } else {
            0.0
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
    fn label(&self, v: &Node) -> String {
        format!("{}", self.emission(v) as char)
    }
}

impl PHMMSampler for LinearPHMM {}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn linear_forward() {
        // generate random seq
        let model = LinearPHMM::from(b"ATCGATTCGATTAGCT");
        let param = crate::hmm::params::PHMMParams::default();
        // let s = b"TTCGGATT";
        let s = b"TTCGATT";
        // let layers = model.forward(&param, s);
        let layers = model.backward(&param, s);

        let kmer_probs: Vec<Vec<Prob>> = layers.iter().map(|layer| layer.to_kmer_prob()).collect();
        for (i, v) in model.nodes.iter().enumerate() {
            print!("{:?}={}", v, model.label(&v));
            for j in 0..kmer_probs.len() {
                print!("\t{}", kmer_probs[j][i].to_log_value());
            }
            print!("\n");
        }
    }
}
