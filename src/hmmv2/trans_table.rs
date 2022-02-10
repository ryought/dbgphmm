//!
//! table structure for storing transition/edge usages
//!
use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use crate::common::Freq;
use crate::prob::Prob;
use crate::vector::{DenseStorage, EdgeVec, Storage};

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
        write!(f, "mm:{}\t", self.mm.to_value())?;
        write!(f, "md:{}\t", self.md.to_value())?;
        write!(f, "im:{}\t", self.im.to_value())?;
        write!(f, "id:{}\t", self.id.to_value())?;
        write!(f, "dm:{}\t", self.dm.to_value())?;
        write!(f, "dd:{}", self.dd.to_value())
    }
}

/// Prob assigned to each edge and each edge types
pub type TransProbs = EdgeVec<DenseStorage<TransProb>>;

/// Frequency (f64) assigned to each edges
pub type EdgeFreqs = EdgeVec<DenseStorage<Freq>>;

//
// visualizer
// TODO move to crate::vector::graph?
//
pub fn draw_edge_vec<N, E, S>(phmm: &PHMMModel<N, E>, ev: &EdgeVec<S>)
where
    N: PHMMNode,
    E: PHMMEdge,
    S: Storage,
    S::Item: std::fmt::Display,
{
    for (e, k, l, _) in phmm.edges() {
        println!(
            "{:?}({}->{})\t{}",
            e,
            phmm.emission(k) as char,
            phmm.emission(l) as char,
            ev[e]
        );
    }
}

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
