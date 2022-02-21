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
        let mut diff = 0f64;
        for (_, d) in self.iter_diff(other) {
            diff += d
        }
        diff
    }
    ///
    /// Calculate `\sum_i |p[i]-p[i]|`
    ///
    pub fn iter_diff<'a, T>(
        &'a self,
        other: &'a Vector<T, Ix>,
    ) -> impl Iterator<Item = (Ix, f64)> + 'a
    where
        T: Storage<Item = Prob>,
    {
        assert_eq!(self.len(), other.len());
        chain!(self.iter().map(|(i, _)| i), other.iter().map(|(i, _)| i))
            .unique()
            .map(move |i| {
                let diff = (self[i].to_value() - other[i].to_value()).abs();
                (i, diff)
            })
    }
    ///
    /// visualize difference of two vectors
    ///
    pub fn show_diff<T>(&self, other: &Vector<T, Ix>)
    where
        T: Storage<Item = Prob>,
    {
        for (i, diff) in self
            .iter_diff(other)
            .sorted_by(|(_, a), (_, b)| b.partial_cmp(a).unwrap())
            .take(100)
        {
            println!("{}: {}", i.index(), diff);
        }
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
        v1[0] = p(0.3);
        let mut v2: Vector<DenseStorage<Prob>, usize> = Vector::new(5, p(0.0));
        v2[1] = p(0.5);
        let d = v1.diff(&v2);
        assert_abs_diff_eq!(d, 0.8);

        let mut v1: Vector<SparseStorage<Prob>, usize> = Vector::new(5, p(0.0));
        v1[0] = p(0.3);
        let mut v2: Vector<SparseStorage<Prob>, usize> = Vector::new(5, p(0.0));
        v2[1] = p(0.5);
        let d = v1.diff(&v2);
        assert_abs_diff_eq!(d, 0.8);

        v1.show_diff(&v2);
    }
}
