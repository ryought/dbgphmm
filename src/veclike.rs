//! Definitions of vec-like structures, that implements `VecLike` trait for setter/gettter
//! to elements
//!
//! ## DenseVec
//! An simple wrapper of std::vec. It will naive example of VecLike implementation.
//!
//! ## SparseVec
//! Vector that
//!
use arrayvec::ArrayVec;
use std::iter::FromIterator;
use std::ops::{Index, IndexMut};

/// VecLike trait abstracts vector-like element access by getter and setter
/// to the index (0 <= index < self.len)
/// with predefined size (it does not support std::vec dynamic sizing)
pub trait VecLike<T: Copy>: Sized {
    fn new(len: usize, value: T) -> Self;
    fn len(&self) -> usize;
    fn get(&self, index: usize) -> T;
    fn set(&mut self, index: usize, value: T);
    fn iter<'a>(&'a self) -> VecLikeIter<'a, T, Self> {
        // VecLike should be Sized in order to pass self here
        VecLikeIter {
            vec: self,
            index: 0,
            phantom: std::marker::PhantomData::<T>,
        }
    }
}

/// Iterator interface
/// This can be created from `VecLike<T>.iter()`
pub struct VecLikeIter<'a, T: Copy, V: VecLike<T>> {
    vec: &'a V,
    index: usize,
    // phantom data has no actual space in VecLikeIter object,
    // but it is necessary for type checking because T should be used somewhere in struct
    phantom: std::marker::PhantomData<T>,
}

impl<'a, T: Copy, V: VecLike<T>> Iterator for VecLikeIter<'a, T, V> {
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.vec.len() {
            let value = self.vec.get(self.index);
            self.index += 1;
            Some(value)
        } else {
            None
        }
    }
}

/// SparseVec max index size parameter
const SIZE: usize = 10;

/// Sparsely storeing vector
/// Space-efficient if there are few non-zero elements in the vector
///
/// For example, `vec[index_0] = value_0` and `vec[index_1] = value_1`
/// SparseVec.index = [index_0, index_1]
/// SparseVec.value = [value_0, value_1]
///
/// When getting non-set index, default_value: T will be returned.
#[derive(Debug, Clone)]
pub struct SparseVec<T: Copy> {
    len: usize,
    default_value: T,
    index: ArrayVec<usize, SIZE>,
    value: ArrayVec<T, SIZE>,
}

/// storing (index, value) in vector for SparseVec
impl<T: Copy> VecLike<T> for SparseVec<T> {
    fn new(len: usize, value: T) -> SparseVec<T> {
        SparseVec {
            len,
            default_value: value,
            index: ArrayVec::<usize, SIZE>::new(),
            value: ArrayVec::<T, SIZE>::new(),
        }
    }
    #[inline(always)]
    fn len(&self) -> usize {
        self.len
    }
    #[inline(always)]
    fn get(&self, index: usize) -> T {
        for i in 0..self.index.len() {
            if self.index[i] == index {
                match self.value.get(i) {
                    Some(value) => return value.clone(),
                    None => panic!(),
                }
            }
        }
        return self.default_value;
    }
    #[inline(always)]
    fn set(&mut self, index: usize, value: T) {
        for i in 0..self.index.len() {
            if self.index[i] == index {
                // update the existing entry
                self.value[i] = value;
                return;
            }
        }

        // add a new entry
        self.index.push(index);
        self.value.push(value);
    }
}

/// Dense vector, a wrapper of std::vec
/// It will use original get, set, len
#[derive(Debug, Clone)]
pub struct DenseVec<T: Copy>(pub Vec<T>);

/// use default std::vec index access for DenseVec
impl<T: Copy> VecLike<T> for DenseVec<T> {
    fn new(len: usize, value: T) -> DenseVec<T> {
        DenseVec(vec![value; len])
    }
    #[inline(always)]
    fn len(&self) -> usize {
        self.0.len()
    }
    #[inline(always)]
    fn get(&self, index: usize) -> T {
        self.0[index]
    }
    #[inline(always)]
    fn set(&mut self, index: usize, value: T) {
        self.0[index] = value
    }
}

impl<T: Copy> FromIterator<T> for DenseVec<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let vec: Vec<T> = iter.into_iter().collect();
        DenseVec(vec)
    }
}

fn head<T: Copy, V: VecLike<T>>(v: &V) -> T {
    v.get(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn veclike_sparse_head() {
        let mut v: SparseVec<usize> = SparseVec::new(10, 0);
        assert_eq!(v.len(), 10);
        assert_eq!(v.get(5), 0);
        v.set(5, 101);
        assert_eq!(v.get(5), 101);
        assert_eq!(head(&v), 0);
        v.set(0, 1111);
        assert_eq!(head(&v), 1111);
    }

    #[test]
    fn veclike_dense_head() {
        let mut v: DenseVec<usize> = DenseVec::new(10, 0);
        assert_eq!(v.len(), 10);
        assert_eq!(v.get(5), 0);
        v.set(5, 101);
        assert_eq!(v.get(5), 101);
        assert_eq!(head(&v), 0);
        v.set(0, 1111);
        assert_eq!(head(&v), 1111);
    }

    #[test]
    fn veclike_sparse_iter() {
        let mut v: SparseVec<usize> = SparseVec::new(10, 0);
        v.set(5, 101);
        v.set(0, 1111);
        v.set(3, 11);
        v.set(7, 89);
        let w: Vec<usize> = v.iter().collect();
        assert_eq!(w, vec![1111, 0, 0, 11, 0, 101, 0, 89, 0, 0]);
    }

    #[test]
    fn veclike_dense_from_iter() {
        let v: Vec<i32> = vec![2, 3, 4, -1, 1, 2, 3];
        let w: DenseVec<i32> = v.into_iter().filter(|&x| x > 0).collect();
        assert_eq!(w.0, vec![2, 3, 4, 1, 2, 3]);
    }
}
