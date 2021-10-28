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
use std::ops::{Index, IndexMut};

/// VecLike trait abstracts vector-like element access by getter and setter
/// to the index (0 <= index < self.size)
/// with predefined size (it does not support std::vec dynamic sizing)
pub trait VecLike<T: Copy> {
    fn new(size: usize, value: T) -> Self;
    fn size(&self) -> usize;
    fn get(&self, index: usize) -> T;
    fn set(&mut self, index: usize, value: T);
}

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
    size: usize,
    default_value: T,
    index: ArrayVec<usize, 10>,
    value: ArrayVec<T, 10>,
}

/// storing (index, value) in vector for SparseVec
impl<T: Copy> VecLike<T> for SparseVec<T> {
    fn new(size: usize, value: T) -> SparseVec<T> {
        SparseVec {
            size,
            default_value: value,
            index: ArrayVec::<usize, 10>::new(),
            value: ArrayVec::<T, 10>::new(),
        }
    }
    fn size(&self) -> usize {
        self.size
    }
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
/// It will use original get, set, size
#[derive(Debug, Clone)]
pub struct DenseVec<T: Copy>(Vec<T>);

/// use default std::vec index access for DenseVec
impl<T: Copy> VecLike<T> for DenseVec<T> {
    fn new(size: usize, value: T) -> DenseVec<T> {
        DenseVec(vec![value; size])
    }
    fn size(&self) -> usize {
        self.0.len()
    }
    fn get(&self, index: usize) -> T {
        self.0[index]
    }
    fn set(&mut self, index: usize, value: T) {
        self.0[index] = value
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
        assert_eq!(v.size(), 10);
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
        assert_eq!(v.size(), 10);
        assert_eq!(v.get(5), 0);
        v.set(5, 101);
        assert_eq!(v.get(5), 101);
        assert_eq!(head(&v), 0);
        v.set(0, 1111);
        assert_eq!(head(&v), 1111);
    }
}
