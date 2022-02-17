//!
//! measureing difference between two vectors
//!
use super::{Indexable, Storage, Vector};
use crate::prob::Prob;
use itertools::{chain, Itertools};

impl<S, Ix> Vector<S, Ix>
where
    S: Storage<Item = Prob>,
    Ix: Indexable,
{
    ///
    /// Log difference score of two vectors whose item is prob.
    ///
    pub fn diff<T>(&self, other: &Vector<T, Ix>) -> f64
    where
        T: Storage<Item = Prob>,
    {
        assert_eq!(self.len(), other.len());
        let mut diff = 0f64;
        for index in chain!(self.iter().map(|(i, _)| i), other.iter().map(|(i, _)| i),).unique() {
            diff += (self[index].to_value() - other[index].to_value()).abs()
        }
        diff
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prob::p;
    use crate::vector::{DenseStorage, SparseStorage};

    #[test]
    fn vector_diff() {
        let mut v1: Vector<DenseStorage<Prob>, usize> = Vector::new(5, p(0.0));
        v1[0] = p(0.5);
        let mut v2: Vector<DenseStorage<Prob>, usize> = Vector::new(5, p(0.0));
        v2[1] = p(0.3);
        let d = v1.diff(&v2);
        assert_abs_diff_eq!(d, 0.8);

        let mut v1: Vector<SparseStorage<Prob>, usize> = Vector::new(5, p(0.0));
        v1[0] = p(0.5);
        let mut v2: Vector<SparseStorage<Prob>, usize> = Vector::new(5, p(0.0));
        v2[1] = p(0.3);
        let d = v1.diff(&v2);
        assert_abs_diff_eq!(d, 0.8);
    }
}
