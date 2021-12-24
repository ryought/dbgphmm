use super::base::{Node, PHMM};
use crate::copy_nums::{Ecdbg, Ncdbg};

impl<'a> PHMM for Ncdbg<'a> {
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
        self.copy_nums.iter().sum::<f64>()
    }
    fn emission(&self, v: &Node) -> u8 {
        self.cdbg.emission(v)
    }
    /// a_vw = c_w / sum of c_w' for w' in childs of v
    fn trans_prob(&self, v: &Node, w: &Node) -> Prob {}
    fn trans_prob(&self, v: &Node, w: &Node) -> Prob {
        match self.cdbg.childs(v).iter().position(|&i| i == *w) {
            Some(index) => self.trans_probs[v.0][index],
            None => Prob::from_prob(0.0),
        }
    }
    fn label(&self, v: &Node) -> String {
        self.cdbg.kmer(v).to_string()
    }
}

impl<'a> PHMM for Ecdbg<'a> {
    /// a_vw = c_vw / c_v
    fn trans_prob(&self, v: &Node, w: &Node) -> Prob {}
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn ncdbg_hmm_0() {}
}
