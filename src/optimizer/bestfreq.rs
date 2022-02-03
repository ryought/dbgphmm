//! Frequency only optimizer, by going up to the best gradient
//!
//! Purpose is same as optimizer::freq, that is
//! (cdbg, freqs) -> (copy_nums)
//! i.e. find the best copy_nums that is closest to freqs
//! by minimizing differences between freqs and copy_nums.

use super::base::ScoreableState;
use super::grad::GDState;
use crate::compressed_dbg::CompressedDBG;
use crate::graph::{IndexedDiGraph, Node};
use std::fmt::Write as FmtWrite;

#[derive(Clone)]
pub struct BestFreqState<'a> {
    cdbg: &'a CompressedDBG,
    idg: &'a IndexedDiGraph,
    freqs: &'a [f64],
    pub copy_nums: Vec<u32>,
}

impl<'a> BestFreqState<'a> {
    pub fn new(
        cdbg: &'a CompressedDBG,
        idg: &'a IndexedDiGraph,
        freqs: &'a [f64],
        copy_nums: Vec<u32>,
    ) -> BestFreqState<'a> {
        BestFreqState {
            cdbg,
            idg,
            freqs,
            copy_nums,
        }
    }
    /// This weight list indicates that how +1/-1 changes of kmer copy number
    /// improves the total difference.
    /// This will be used as a edge weight in min-mean-weight-cycle problem.
    /// Each kmer has two entries:
    /// - `weights[2*i] = (score gain of +1 in i-th kmer)`
    /// - `weights[2*i + 1] = (score gain of -1 in i-th kmer)`
    fn weights(&self) -> Vec<f64> {
        (0..self.cdbg.n_kmers())
            .flat_map(|i| {
                let freq = self.freqs[i];
                let copy_num = self.copy_nums[i];

                // current score
                let dist_now = (freq - copy_num as f64).powi(2);

                // next score if increased or decreased
                let dist_inc = (freq - (copy_num as f64 + 1.0)).powi(2);
                let dist_dec = (freq - (copy_num as f64 - 1.0)).powi(2);

                if !self.cdbg.is_emitable(&Node(i)) {
                    // do not care score of non-emittable nodes
                    if copy_num > 0 {
                        return vec![0.0, 0.0];
                    } else {
                        return vec![0.0, f64::INFINITY];
                    }
                }

                // if copy_num == 0, the kmer is not down-movable
                if copy_num > 0 {
                    vec![dist_inc - dist_now, dist_dec - dist_now]
                } else {
                    vec![dist_inc - dist_now, f64::INFINITY]
                }
            })
            .collect()
    }

    /// Create a new copy_nums
    /// by updating along the min-mean-weight-cycle
    fn best_updated_copy_nums(&self) -> Option<Vec<u32>> {
        let weights = self.weights();
        match self.idg.minimum_mean_weight_cycle(&Node(0), &weights) {
            Some(cycle) => {
                let mut copy_nums = self.copy_nums.clone();
                for e in cycle.iter() {
                    if e.0 % 2 == 0 {
                        // +1
                        copy_nums[e.0 / 2] += 1
                    } else {
                        // -1
                        copy_nums[e.0 / 2] -= 1
                    }
                }
                Some(copy_nums)
            }
            None => None,
        }
    }
}

impl<'a> ScoreableState for BestFreqState<'a> {
    /// calc score of the state
    fn score(&self) -> f64 {
        let diff: f64 = self
            .copy_nums
            .iter()
            .zip(self.freqs.iter())
            .enumerate()
            .filter(|&(i, _)| self.cdbg.is_emitable(&Node(i)))
            .map(|(_, x)| x)
            .map(|(&cn, &f)| (cn as f64 - f).powi(2))
            .sum();
        diff * -1f64
    }
    fn fill_score(&mut self) -> f64 {
        // not using a cache, so simply returns score()
        self.score()
    }
    #[allow(unused_must_use)]
    fn as_string(&self) -> String {
        let mut s = String::new();
        let freq_strs: Vec<_> = self
            .freqs
            .iter()
            .map(|freq| format!("{:.2}", freq))
            .collect();
        write!(
            &mut s,
            "{}\t{}\t{}\t{}\n{:?}\t[{}]",
            self.score(),
            0.0,
            self.cdbg.total_emitable_copy_num(&self.copy_nums),
            self.cdbg.to_seqs_string(&self.copy_nums),
            self.copy_nums,
            freq_strs.join(", "),
        );
        s
    }
}

impl<'a> GDState for BestFreqState<'a> {
    /// calc the current difference between freq and copy_nums
    /// and find the best cycle that improves the diffs
    /// by min-mean-weight-cycle.
    /// The number of neighbors is always 0-or-1.
    fn neighbors(&self) -> Vec<BestFreqState<'a>> {
        match self.best_updated_copy_nums() {
            // if it can be updated, create new BestFreqState
            Some(copy_nums) => {
                vec![BestFreqState {
                    cdbg: self.cdbg,
                    idg: self.idg,
                    freqs: self.freqs,
                    copy_nums,
                }]
            }
            // if no updates, return empty vector
            None => Vec::new(),
        }
    }
    fn is_duplicate(&self, other: &BestFreqState) -> bool {
        self.copy_nums == other.copy_nums
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::optimizer::grad::GradientDescent;

    #[test]
    fn simple_best_freq_state() {
        let seqs = vec![b"ATTCGATCGATTT".to_vec()];
        let (cdbg, copy_nums) = CompressedDBG::from_seqs(&seqs, 8);
        let freqs = cdbg.copy_nums_to_freqs(&copy_nums);
        let idg = cdbg.to_indexed_digraph();
        // let s = BestFreqState::new(&cdbg, &idg, &freqs, copy_nums);

        // start from [0,0,0,....]
        let s = BestFreqState::new(&cdbg, &idg, &freqs, vec![0; cdbg.n_kmers()]);
        assert_eq!(s.score(), -20.0);
        for (i, w) in s.weights().into_iter().enumerate() {
            if i % 2 == 0 {
                assert_eq!(w, -1.0);
            } else {
                assert_eq!(w, f64::INFINITY);
            }
        }
        let copy_nums_new = s.best_updated_copy_nums();
        assert_eq!(copy_nums_new, Some(vec![1; cdbg.n_kmers()]));

        let neighbors = s.neighbors();
        assert_eq!(neighbors.len(), 1);
        assert_eq!(neighbors[0].score(), 0.0);
    }

    #[test]
    fn simple_best_freq_grad() {
        let seqs = vec![b"ATTCGATCGATTT".to_vec()];
        let (cdbg, copy_nums) = CompressedDBG::from_seqs(&seqs, 8);
        let freqs = cdbg.copy_nums_to_freqs(&copy_nums);
        let idg = cdbg.to_indexed_digraph();
        let s = BestFreqState::new(&cdbg, &idg, &freqs, vec![0; cdbg.n_kmers()]);

        let g = GradientDescent::new(100, true);
        g.run(s);
    }
}
