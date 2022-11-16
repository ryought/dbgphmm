//!
//! Greedy search of posterior probability
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase, EdgeCopyNums, NodeCopyNums};
use crate::common::CopyNum;
use crate::dbg::cycle::CopyNumsUpdateInfo;
use crate::e2e::Dataset;
use crate::greedy::{GreedyInstance, GreedyScore, GreedySearcher};
use crate::hmmv2::params::PHMMParams;
use crate::kmer::kmer::KmerLike;
use crate::prob::Prob;

//
// Instance
//
#[derive(Clone)]
pub struct DbgCopyNumsInstance<K: KmerLike> {
    copy_nums: NodeCopyNums,
    info: CopyNumsUpdateInfo<K>,
}
impl<K: KmerLike> DbgCopyNumsInstance<K> {
    fn new(copy_nums: NodeCopyNums, info: CopyNumsUpdateInfo<K>) -> Self {
        DbgCopyNumsInstance { copy_nums, info }
    }
    fn copy_nums(&self) -> &NodeCopyNums {
        &self.copy_nums
    }
    fn info(&self) -> &CopyNumsUpdateInfo<K> {
        &self.info
    }
}
impl<K: KmerLike> GreedyInstance for DbgCopyNumsInstance<K> {
    type Key = NodeCopyNums;
    fn key(&self) -> &NodeCopyNums {
        &self.copy_nums
    }
}

//
// Score
//
#[derive(Clone, Copy)]
pub struct DbgCopyNumsScore {
    /// Likelihood P(R|G)
    pub p_rg: Prob,
    /// Prior P(G)
    pub p_g: Prob,
}
impl DbgCopyNumsScore {
    pub fn new(p_rg: Prob, p_g: Prob) -> Self {
        DbgCopyNumsScore { p_rg, p_g }
    }
}
impl GreedyScore for DbgCopyNumsScore {
    /// Unnormalized posterior P(R,G) = P(R|G)P(G)
    fn prob(&self) -> Prob {
        self.p_rg * self.p_g
    }
}

pub type Posterior<K> = Vec<(Prob, DbgCopyNumsInstance<K>, DbgCopyNumsScore)>;

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    pub fn search_posterior(
        &self,
        dataset: &Dataset,
        max_neighbor_depth: usize,
        max_move: usize,
        genome_size_expected: CopyNum,
        genome_size_sigma: CopyNum,
    ) -> Posterior<N::Kmer> {
        let instance_init =
            DbgCopyNumsInstance::new(self.to_node_copy_nums(), CopyNumsUpdateInfo::empty());
        let mut searcher = GreedySearcher::new(
            instance_init,
            |instance| {
                eprintln!("[to_score] calculating prob {}", instance.copy_nums());
                let mut dbg = self.clone();
                dbg.set_node_copy_nums(instance.copy_nums());
                let p_rg = dbg.to_full_prob(dataset.params(), dataset.reads());
                let p_g = dbg.to_prior_prob(genome_size_expected, genome_size_sigma);
                eprintln!("[to_score] P(R|G)={} P(G)={}", p_rg, p_g);
                DbgCopyNumsScore::new(p_rg, p_g)
            },
            |instance| {
                eprintln!("[to_neighbors] calculating neighbors...");
                let mut dbg = self.clone();
                dbg.set_node_copy_nums(instance.copy_nums());
                let neighbors: Vec<_> = dbg
                    .neighbor_copy_nums_fast_compact_with_info(max_neighbor_depth)
                    .into_iter()
                    .filter(|(copy_nums, _)| copy_nums.sum() > 0) // remove null genome
                    .map(|(copy_nums, info)| DbgCopyNumsInstance::new(copy_nums, info))
                    .collect();
                eprintln!("[to_neighbors] #neighbors={}", neighbors.len());
                neighbors
            },
        );

        searcher.search(max_move);
        searcher.into_posterior_distribution()
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::e2e::{generate_simple_genome_mock, generate_small_tandem_repeat};

    #[test]
    fn greedy_simple_genome() {
        let experiment = generate_simple_genome_mock();
        let dbg_draft = experiment.dbg_draft.clone().unwrap();
        let copy_nums = dbg_draft.neighbor_copy_nums_fast_compact(100);

        for copy_num in copy_nums.iter() {
            println!("c={}", copy_num);
        }

        let sigma = 10;
        let s = dbg_draft.search_posterior(
            experiment.dataset(),
            100,
            10,
            experiment.genome_size(),
            sigma,
        );
        for (p_gr, instance, score) in s.into_iter() {
            println!(
                "P(G|R)={} (P(R|G)={}, P(G)={}) {} {}",
                p_gr,
                score.p_rg,
                score.p_g,
                instance.copy_nums(),
                instance.info(),
            );
        }
    }
    #[test]
    fn greedy_tandem_repeat() {
        let experiment = generate_small_tandem_repeat();
        let dbg_draft_true = experiment.dbg_draft_true.clone().unwrap();
        let copy_nums = dbg_draft_true.neighbor_copy_nums_fast_compact(5);

        println!("n_neighbors={}", copy_nums.len());

        let approx = experiment.dbg_draft.clone().unwrap().to_node_copy_nums();
        let copy_nums_true = dbg_draft_true.to_node_copy_nums();
        println!("diff(approx, true)={}", approx.dist(&copy_nums_true));

        let sigma = 100;
        let distribution = dbg_draft_true.search_posterior(
            experiment.dataset(),
            5,
            1,
            experiment.genome_size(),
            sigma,
        );

        for (p_gr, instance, score) in distribution.iter() {
            println!(
                "P(G|R)={} (P(R|G)={}, P(G)={}) {} {} {}",
                p_gr,
                score.p_rg,
                score.p_g,
                instance.info(),
                instance.copy_nums().dist(&copy_nums_true),
                instance.copy_nums(),
            );
        }
        let neighbors: Vec<_> = distribution
            .iter()
            .map(|(p_gr, instance, _score)| (instance.copy_nums().clone(), *p_gr))
            .collect();
        dbg_draft_true.to_kmer_distribution(&neighbors);
    }
}
