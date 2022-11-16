//!
//! Greedy searcher
//!
//! * Init: start from initial state
//! * Neighbor: calculate scores of neighboring states.
//! * Move: pick a state with the highest score / pick a unexamined state satisfying certain conditions
//! * Stop: there is no state with higher score / there is no move
//!
use crate::prob::Prob;
use fnv::FnvHashMap as HashMap;
use itertools::Itertools;
use std::hash::Hash;

///
///
///
pub trait GreedyInstance: Clone {
    type Key: Clone + Eq + Hash;
    ///
    /// a unique key of the instance
    ///
    fn key(&self) -> &Self::Key;
    // ///
    // /// generate neighboring instances
    // ///
    // fn neighbors(&self) -> Vec<Self>;
    // /// calculate unnormalized posterior probability `P(R|G)P(G)`
    // ///
    // /// TODO
    // /// * we want to pass a reference to data (which is necessary for calculating the score of the
    // /// state).
    // fn score(&self) -> Prob;
}

pub trait GreedyScore: Copy {
    fn prob(&self) -> Prob;
}

///
/// * pick a highest score state efficiently
///
pub struct GreedySearcher<I: GreedyInstance, S: GreedyScore, F: Fn(&I) -> S, G: Fn(&I) -> Vec<I>> {
    // TODO to_score and to_neighbors should be passed as args?
    to_score: F,
    to_neighbors: G,
    current_instance: I,
    ///
    /// Sum of probabilities of items stored in `self.instances`
    ///
    total_prob: Prob,
    instances: HashMap<I::Key, (I, S)>,
    // history: Vec<(I, Prob)>,
}

impl<I: GreedyInstance, S: GreedyScore, F: Fn(&I) -> S, G: Fn(&I) -> Vec<I>>
    GreedySearcher<I, S, F, G>
{
    ///
    /// create GreedySearch from a start point initial instance.
    ///
    pub fn new(instance_init: I, to_score: F, to_neighbors: G) -> Self {
        let p = to_score(&instance_init);
        let mut instances = HashMap::default();
        instances.insert(instance_init.key().clone(), (instance_init.clone(), p));
        GreedySearcher {
            to_score,
            to_neighbors,
            current_instance: instance_init,
            total_prob: p.prob(),
            instances,
        }
    }
    ///
    /// calculate score of an instance
    ///
    fn calc_score(&self, instance: &I) -> S {
        (self.to_score)(instance)
    }
    ///
    /// calculate neighbors of an instance
    ///
    fn calc_neighbors(&self, instance: &I) -> Vec<I> {
        (self.to_neighbors)(instance)
    }
    ///
    ///
    fn get_highest_instance(&self) -> (&I, S) {
        self.instances
            .values()
            .max_by_key(|(_, p)| p.prob())
            .map(|(instance, p)| (instance, *p))
            .unwrap()
    }
    ///
    ///
    ///
    pub fn search_once(&mut self) -> usize {
        // find neighbor
        let neighbors = self.calc_neighbors(&self.current_instance);
        let mut n_new_neighbors = 0;
        // calculate score of new neighbors
        for neighbor in neighbors.into_iter() {
            if !self.instances.contains_key(neighbor.key()) {
                let p = self.calc_score(&neighbor);
                self.instances.insert(neighbor.key().clone(), (neighbor, p));
                self.total_prob += p.prob();
                n_new_neighbors += 1;
            }
        }
        // store the highest score element in current
        let (next_instance, _) = self.get_highest_instance();
        self.current_instance = next_instance.clone();

        n_new_neighbors
    }
    ///
    /// run `search_once` recursively
    ///
    pub fn search(&mut self, max_move: usize) {
        for i in 0..max_move {
            let n_new_neighbors = self.search_once();
            if n_new_neighbors == 0 {
                break;
            }
        }
    }
    ///
    ///
    ///
    pub fn into_posterior_distribution(self) -> Vec<(Prob, I, S)> {
        // z: normalizing constant
        let z = self.total_prob;
        self.instances
            .into_iter()
            .map(|(_, (instance, score))| (score.prob() / z, instance, score))
            .sorted_by_key(|(p, _, _)| *p)
            .collect()
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aa() {}
}
