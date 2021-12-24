use super::base::{Node, PHMM};
use crate::copy_nums::{Ecdbg, Ncdbg};
use crate::prob::Prob;

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
        self.cdbg.is_adjacent(v, w)
    }
    fn copy_num(&self, v: &Node) -> f64 {
        self.copy_num(v) as f64
    }
    fn total_copy_num(&self) -> f64 {
        let total_copy_num = self.cdbg.total_emitable_copy_num(&self.copy_nums.0);
        total_copy_num as f64
    }
    fn emission(&self, v: &Node) -> u8 {
        self.cdbg.emission(v)
    }
    /// a_vw = c_w / sum of c_w' for w' in childs of v
    fn trans_prob(&self, v: &Node, w: &Node) -> Prob {
        let total_cn_of_childs: u32 = self
            .cdbg
            .childs(&v)
            .iter()
            .map(|&w| {
                if self.cdbg.is_emitable(&w) {
                    self.copy_num(&w)
                } else {
                    0
                }
            })
            .sum();
        let cn_of_child = self.copy_num(w);

        if self.is_emitable(&v) && total_cn_of_childs > 0 {
            Prob::from_prob(f64::from(cn_of_child) / f64::from(total_cn_of_childs))
        } else {
            Prob::from_prob(0.0)
        }
    }
    fn label(&self, v: &Node) -> String {
        self.cdbg.kmer(v).to_string()
    }
}

impl<'a> PHMM for Ecdbg<'a> {
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
        self.cdbg.is_adjacent(v, w)
    }
    fn copy_num(&self, v: &Node) -> f64 {
        self.copy_num(v) as f64
    }
    fn total_copy_num(&self) -> f64 {
        let total_copy_num = self.cdbg.total_emitable_copy_num(&self.node_copy_nums.0);
        total_copy_num as f64
    }
    fn emission(&self, v: &Node) -> u8 {
        self.cdbg.emission(v)
    }
    /// a_vw = c_vw / c_v
    fn trans_prob(&self, v: &Node, w: &Node) -> Prob {
        let copy_num_v = self.copy_num(v);
        let copy_num_vw = self.edge_copy_num(v, w);

        if self.is_emitable(&v) && copy_num_v > 0 {
            Prob::from_prob(f64::from(copy_num_vw) / f64::from(copy_num_v))
        } else {
            Prob::from_prob(0.0)
        }
    }
    fn label(&self, v: &Node) -> String {
        self.cdbg.kmer(v).to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mocks::test_cdbg_01;
    #[test]
    fn ncdbg_hmm_0() {
        let (cdbg, cn, _) = test_cdbg_01();
        let ncdbg = Ncdbg::new(&cdbg, cn);
        println!("{}", ncdbg.as_dot());
    }
    #[test]
    fn ecdbg_hmm_0() {
        let (cdbg, cn, ecn) = test_cdbg_01();
        let ecdbg = Ecdbg::new(&cdbg, cn, ecn);
        println!("{}", ecdbg.as_dot());
    }
}
