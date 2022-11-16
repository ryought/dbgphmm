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
use std::hash::Hash;

///
///
///
pub trait GreedyInstance: Clone + Eq + Hash {
    // ///
    // /// generate neighboring instances
    // ///
    // fn neighbors(&self) -> Vec<Self>;
    // type Key: Clone;
    // ///
    // /// a unique key of the instance
    // ///
    // fn key(&self) -> &Self::Key;
    // /// calculate unnormalized posterior probability `P(R|G)P(G)`
    // ///
    // /// TODO
    // /// * we want to pass a reference to data (which is necessary for calculating the score of the
    // /// state).
    // fn score(&self) -> Prob;
}

///
/// * pick a highest score state efficiently
///
pub struct GreedySearcher<I: GreedyInstance, F: Fn(&I) -> Prob, G: Fn(&I) -> Vec<I>> {
    to_score: F,
    to_neighbors: G,
    current_instance: I,
    // history: Vec<(I, Prob)>,
    total_prob: Prob,
    instances: HashMap<I, Prob>,
}

impl<I: GreedyInstance, F: Fn(&I) -> Prob, G: Fn(&I) -> Vec<I>> GreedySearcher<I, F, G> {
    ///
    /// create GreedySearch from a start point initial instance.
    ///
    pub fn new(instance_init: I, to_score: F, to_neighbors: G) -> Self {
        let p = to_score(&instance_init);
        let mut instances = HashMap::default();
        instances.insert(instance_init.clone(), p);
        GreedySearcher {
            to_score,
            to_neighbors,
            current_instance: instance_init.clone(),
            total_prob: p,
            instances,
        }
    }
    ///
    /// calculate score of an instance
    ///
    fn calc_score(&self, instance: &I) -> Prob {
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
    fn get_highest_instance(&self) -> (&I, Prob) {
        self.instances
            .iter()
            .max_by_key(|(_, p)| *p)
            .map(|(instance, &p)| (instance, p))
            .unwrap()
    }
    ///
    ///
    pub fn search_once(&mut self) -> usize {
        // find neighbor
        let neighbors = self.calc_neighbors(&self.current_instance);
        let mut n_new_neighbors = 0;
        // calculate score of new neighbors
        for neighbor in neighbors.into_iter() {
            if !self.instances.contains_key(&neighbor) {
                let p = self.calc_score(&neighbor);
                self.instances.insert(neighbor, p);
                self.total_prob += p;
                n_new_neighbors += 1;
            }
        }
        // store the highest score element in current
        let (next_instance, _) = self.get_highest_instance();
        self.current_instance = next_instance.clone();

        n_new_neighbors
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
