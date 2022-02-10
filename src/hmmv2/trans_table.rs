//!
//! table structure for storing transition/edge usages
//!
use crate::common::Freq;
use crate::prob::Prob;
use crate::vector::{DenseStorage, EdgeVec, NodeVec, Storage};

#[derive(Debug, Copy, PartialEq, Clone)]
pub struct TransProb {
    /// Mk -> Ml transition
    pub mm: Prob,
    /// Mk -> Dl transition
    pub md: Prob,
    /// Ik -> Ml transition
    pub im: Prob,
    /// Ik -> Dl transition
    pub id: Prob,
    /// Dk -> Ml transition
    pub dm: Prob,
    /// Dk -> Dl transition
    pub dd: Prob,
}

impl TransProb {
    /// Constructor of TransProb with all p=0
    pub fn zero() -> Self {
        TransProb {
            mm: Prob::from_prob(0.0),
            md: Prob::from_prob(0.0),
            im: Prob::from_prob(0.0),
            id: Prob::from_prob(0.0),
            dm: Prob::from_prob(0.0),
            dd: Prob::from_prob(0.0),
        }
    }
    /// sum of probabilities of all elements
    /// That is `mm + md + im + id + dm + dd`
    pub fn sum(&self) -> Prob {
        self.mm + self.md + self.im + self.id + self.dm + self.dd
    }
}

impl std::fmt::Display for TransProb {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "m -> m\t{}", self.mm)?;
        writeln!(f, "m -> d\t{}", self.md)?;
        writeln!(f, "i -> m\t{}", self.im)?;
        writeln!(f, "i -> d\t{}", self.id)?;
        writeln!(f, "d -> m\t{}", self.dm)?;
        writeln!(f, "d -> d\t{}", self.dd)
    }
}

/// Prob assigned to each edge and each edge types
pub type TransProbs = EdgeVec<DenseStorage<TransProb>>;

/// Frequency (f64) assigned to each edges
pub type EdgeFreqs = EdgeVec<DenseStorage<Freq>>;

// test
//
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trans_prob_sum() {
        let mut tp = TransProb::zero();
        println!("{}", tp);
        assert_eq!(tp.sum(), Prob::from_prob(0.0));
        tp.mm = Prob::from_prob(0.2);
        tp.md = Prob::from_prob(0.1);
        tp.im = Prob::from_prob(0.3);
        tp.id = Prob::from_prob(0.4);
        tp.dm = Prob::from_prob(0.2);
        tp.dd = Prob::from_prob(0.01);
        assert_abs_diff_eq!(
            tp.sum(),
            Prob::from_prob(0.2 + 0.1 + 0.3 + 0.4 + 0.2 + 0.01)
        );
        println!("{}", tp);
    }
}
