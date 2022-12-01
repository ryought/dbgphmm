//!
//! Greedy search of posterior probability
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase, EdgeCopyNums, NodeCopyNums};
use crate::common::CopyNum;
use crate::dbg::neighbor::CopyNumsUpdateInfo;
use crate::dbg::phmm::EvalResult;
use crate::e2e::Dataset;
use crate::greedy::{GreedyInstance, GreedyScore, GreedySearcher};
use crate::hmmv2::params::PHMMParams;
use crate::kmer::kmer::KmerLike;
use crate::prob::Prob;
use fnv::FnvHashSet as HashSet;
use itertools::Itertools;
use std::time::{Duration, Instant};

//
// Instance
//
#[derive(Clone)]
pub struct DbgCopyNumsInstance<K: KmerLike> {
    copy_nums: NodeCopyNums,
    info: CopyNumsUpdateInfo<K>,
    move_count: usize,
}
impl<K: KmerLike> DbgCopyNumsInstance<K> {
    pub fn new(copy_nums: NodeCopyNums, info: CopyNumsUpdateInfo<K>, move_count: usize) -> Self {
        DbgCopyNumsInstance {
            copy_nums,
            info,
            move_count,
        }
    }
    pub fn copy_nums(&self) -> &NodeCopyNums {
        &self.copy_nums
    }
    pub fn info(&self) -> &CopyNumsUpdateInfo<K> {
        &self.info
    }
    pub fn move_count(&self) -> usize {
        self.move_count
    }
}
impl<K: KmerLike> GreedyInstance for DbgCopyNumsInstance<K> {
    type Key = NodeCopyNums;
    fn key(&self) -> &NodeCopyNums {
        &self.copy_nums
    }
}

impl GreedyScore for EvalResult {
    fn prob(&self) -> Prob {
        self.posterior()
    }
}

pub type Posterior<K> = Vec<(Prob, DbgCopyNumsInstance<K>, EvalResult)>;

pub fn get_max_posterior_instance<K: KmerLike>(
    posterior: &Posterior<K>,
) -> &DbgCopyNumsInstance<K> {
    posterior
        .iter()
        .max_by_key(|(p, _, _)| *p)
        .map(|(_, instance, _)| instance)
        .unwrap()
}

//
// Posterior
//
#[derive(Clone)]
pub struct DbgCopyNumsPosterior<N: DbgNode, E: DbgEdge> {
    dbg: Dbg<N, E>,
    copy_nums: Vec<(Prob, DbgCopyNumsInstance<N::Kmer>, EvalResult)>,
}
impl<N: DbgNode, E: DbgEdge> DbgCopyNumsPosterior<N, E> {}

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    /// search copy_nums with extending k
    ///
    pub fn search_posterior_with_extension(
        &self,
        dataset: &Dataset,
        max_neighbor_depth: usize,
        max_move: usize,
        genome_size_expected: CopyNum,
        genome_size_sigma: CopyNum,
    ) -> Posterior<N::Kmer> {
        unimplemented!();
    }
    ///
    ///
    ///
    pub fn search_posterior<F>(
        &self,
        dataset: &Dataset,
        max_neighbor_depth: usize,
        max_move: usize,
        genome_size_expected: CopyNum,
        genome_size_sigma: CopyNum,
        on_move: F,
    ) -> Posterior<N::Kmer>
    where
        F: Fn(&DbgCopyNumsInstance<N::Kmer>),
    {
        let instance_init =
            DbgCopyNumsInstance::new(self.to_node_copy_nums(), CopyNumsUpdateInfo::empty(), 0);
        eprintln!("creating hint");
        let reads_with_hints = self.generate_hints(dataset.reads(), dataset.params());
        let mut searcher = GreedySearcher::new(
            instance_init,
            |instance| {
                let mut dbg = self.clone();
                dbg.set_node_copy_nums(instance.copy_nums());
                let r = dbg.evaluate_with_hint(
                    dataset.params(),
                    &reads_with_hints,
                    genome_size_expected,
                    genome_size_sigma,
                );
                // eprintln!("# [to_score/#{}] {}", instance.move_count(), r);
                r
            },
            |instance| {
                on_move(instance);
                let mut dbg = self.clone();
                let start = Instant::now();
                dbg.set_node_copy_nums(instance.copy_nums());
                let neighbors: Vec<_> = dbg
                    .neighbor_copy_nums_fast_compact_with_info(max_neighbor_depth, false)
                    .into_iter()
                    .filter(|(copy_nums, _)| copy_nums.sum() > 0) // remove null genome
                    .map(|(copy_nums, info)| {
                        DbgCopyNumsInstance::new(copy_nums, info, instance.move_count + 1)
                    })
                    .collect();
                let duration = start.elapsed();
                // eprintln!(
                //     "[to_neighbors/#{}] found {} neighbors (in {} ms)",
                //     instance.move_count(),
                //     neighbors.len(),
                //     duration.as_millis(),
                // );
                neighbors
            },
        );

        searcher.search(max_move);
        searcher.into_posterior_distribution()
    }
    ///
    ///
    ///
    pub fn search_posterior_with_restriction(
        &self,
        dataset: &Dataset,
        max_neighbor_depth: usize,
        max_move: usize,
        genome_size_expected: CopyNum,
        genome_size_sigma: CopyNum,
    ) -> Posterior<N::Kmer> {
        let instance_init =
            DbgCopyNumsInstance::new(self.to_node_copy_nums(), CopyNumsUpdateInfo::empty(), 0);
        let mut searcher = GreedySearcher::new(
            instance_init,
            |instance| {
                let mut dbg = self.clone();
                dbg.set_node_copy_nums(instance.copy_nums());
                let r = dbg.evaluate(
                    dataset.params(),
                    dataset.reads(),
                    genome_size_expected,
                    genome_size_sigma,
                );
                eprintln!("[to_score/#{}] {}", instance.move_count(), r);
                r
            },
            |instance| {
                let mut dbg = self.clone();
                let start = Instant::now();
                dbg.set_node_copy_nums(instance.copy_nums());
                let neighbors_raw: Vec<_> = dbg
                    .neighbor_copy_nums_fast_compact_with_info(max_neighbor_depth, true)
                    .into_iter()
                    .filter(|(copy_nums, _)| copy_nums.sum() > 0) // remove null genome
                    .collect();

                let mut zero_one_neighbors = HashSet::default();
                let mut neighbors = Vec::new();
                for (neighbor, info) in neighbors_raw.into_iter() {
                    let zero_one = neighbor.zero_one();
                    if !zero_one_neighbors.contains(&zero_one) {
                        zero_one_neighbors.insert(zero_one);
                        neighbors.push((neighbor, info));
                    }
                }

                let neighbors: Vec<_> = neighbors
                    .into_iter()
                    .map(|(copy_nums, info)| {
                        DbgCopyNumsInstance::new(copy_nums, info, instance.move_count + 1)
                    })
                    .collect();
                let duration = start.elapsed();
                eprintln!(
                    "[to_neighbors/#{}] found {} neighbors (in {} ms)",
                    instance.move_count(),
                    neighbors.len(),
                    duration.as_millis(),
                );
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
    use crate::dbg::SimpleDbg;
    use crate::e2e::{
        generate_difficult_diploid_tandem_repeat_dataset, generate_simple_genome_mock,
        generate_small_tandem_repeat,
    };
    use crate::kmer::VecKmer;

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
            |_| {},
        );
        for (p_gr, instance, score) in s.into_iter() {
            println!(
                "P(G|R)={} (P(R|G)={}, P(G)={}) {} {}",
                p_gr,
                score.p_rg(),
                score.p_g(),
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
            |_| {},
        );

        for (p_gr, instance, score) in distribution.iter() {
            println!(
                "P(G|R)={} (P(R|G)={}, P(G)={}) {} {} {} {}",
                p_gr,
                score.p_rg(),
                score.p_g(),
                instance.move_count(),
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
    #[test]
    #[ignore = "takes long time (~1 hour)"]
    fn dbg_sample_posterior_for_difficult_tandem_repeat() {
        let dataset = generate_difficult_diploid_tandem_repeat_dataset();
        let mut dbg: SimpleDbg<VecKmer> =
            SimpleDbg::create_draft_from_fragment_seqs_with_adjusted_coverage(
                12,
                dataset.reads(),
                dataset.coverage(),
                dataset.reads().average_length(),
                dataset.params().p_error().to_value(),
            );
        let (copy_nums_true, _) = dbg.to_copy_nums_of_styled_seqs(dataset.genome()).unwrap();
        let copy_nums_draft = dbg.to_node_copy_nums();
        let distribution =
            dbg.search_posterior_with_restriction(&dataset, 5, 10, dataset.genome_size(), 100);

        for (p_gr, instance, score) in distribution.iter().sorted_by_key(|(p, _, _)| *p) {
            println!(
                "P(G|R)={} (P(R|G)={}, P(G)={}) {} {}",
                p_gr,
                score.p_rg(),
                score.p_g(),
                instance.move_count(),
                instance.copy_nums().dist(&copy_nums_true),
            );
        }
    }
}
