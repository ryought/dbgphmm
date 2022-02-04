//!
//! Sparse storage that uses `ArrayVec`
//!
use super::{IterableStorage, Storage, Vector};
use arrayvec::ArrayVec;

/// SparseStorage max index size parameter
const SIZE: usize = 200;

/// Sparse storage powered by `ArrayVec`
/// A type of items is represented as T.
#[derive(Debug, Clone)]
pub struct SparseStorage<T> {
    /// virtual size of this storage
    size: usize,
    /// default value
    /// The value of unused index wiil be filled with
    /// this `default_value`.
    default_value: T,
    /// ArrayVec of elements `(index, value)`
    elements: ArrayVec<(usize, T), SIZE>,
}

impl<T> Storage for SparseStorage<T>
where
    T: Copy,
{
    type Item = T;
    fn new(size: usize, default_value: T) -> SparseStorage<T> {
        SparseStorage {
            size,
            default_value,
            elements: ArrayVec::<(usize, T), SIZE>::new(),
        }
    }
    #[inline(always)]
    fn size(&self) -> usize {
        self.size
    }
    fn get(&self, index: usize) -> &T {
        assert!(index < self.size());

        // search for existing entries
        for i in 0..self.elements.len() {
            if self.elements[i].0 == index {
                return &self.elements[i].1;
            }
        }

        // not found
        &self.default_value
    }
    fn get_mut(&mut self, index: usize) -> &mut T {
        assert!(index < self.size());

        // search for existing entries
        for i in 0..self.elements.len() {
            if self.elements[i].0 == index {
                return &mut self.elements[i].1;
            }
        }

        // add a new entry and return the reference to it
        self.elements.push((index, self.default_value));
        let n = self.elements.len();
        return &mut self.elements[n - 1].1;
    }
}

impl<'a, T> IterableStorage<'a> for SparseStorage<T>
where
    T: Copy + 'a,
{
    type IndexIterator = SparseStorageIterator<'a, T>;
    fn indexiter(&'a self) -> Self::IndexIterator {
        SparseStorageIterator {
            i: 0,
            storage: &self,
        }
    }
}

/// Iterator on SparseStorage (index, item)
pub struct SparseStorageIterator<'a, T: Copy + 'a> {
    /// current index on the arrayvec on the sparse storage
    i: usize,
    /// reference to the storage
    storage: &'a SparseStorage<T>,
}

impl<'a, T: Copy + 'a> Iterator for SparseStorageIterator<'a, T> {
    type Item = (usize, T);
    fn next(&mut self) -> Option<Self::Item> {
        if self.i < self.storage.elements.len() {
            let i = self.i;
            let (index, value) = self.storage.elements[i];
            self.i += 1;
            Some((index, value))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sparse_storage() {
        // u32
        let mut s: SparseStorage<u32> = SparseStorage::new(10, 0);
        *s.get_mut(5) = 111;
        *s.get_mut(3) = 22;
        assert_eq!(*s.get(0), 0);
        assert_eq!(*s.get(3), 22);
        assert_eq!(*s.get(5), 111);
        assert_eq!(s.size(), 10);

        // clone
        let mut s2 = s.clone();
        assert_eq!(*s2.get(3), 22);
        *s2.get_mut(3) = 21;
        assert_eq!(*s2.get(3), 21);
        assert_eq!(*s.get(3), 22);

        // f64
        let mut s: SparseStorage<f64> = SparseStorage::new(10, 0.0);
        *s.get_mut(5) = 12.11;
        *s.get_mut(3) = 10.0;
        assert_eq!(*s.get(0), 0.0);
        assert_eq!(*s.get(3), 10.0);
        assert_eq!(*s.get(5), 12.11);
    }
    #[test]
    #[should_panic]
    fn sparse_storage_outside() {
        let mut s: SparseStorage<u32> = SparseStorage::new(3, 0);
        *s.get_mut(3) = 22;
    }
    #[test]
    fn sparse_storage_iter() {
        let mut s: SparseStorage<u32> = SparseStorage::new(4, 0);
        *s.get_mut(0) = 111;
        *s.get_mut(2) = 10;
        let v: Vec<(usize, u32)> = s.indexiter().collect();
        assert_eq!(v, vec![(0, 111), (2, 10)]);
    }
    #[test]
    fn sparse_storage_vector() {
        let mut v: Vector<SparseStorage<u32>> = Vector::new(5, 0);
        v[0] = 100;
        v[3] = 222;
        assert_eq!(v[0], 100);
        assert_eq!(v[1], 0);
        let w: Vec<(usize, u32)> = v.iter().collect();
        println!("{:?}", v);
        println!("{:?}", w);
        assert_eq!(w, vec![(0, 100), (3, 222)]);
    }
    #[test]
    fn sparse_storage_vector_add_mul() {
        let mut v1: Vector<SparseStorage<u32>> = Vector::new(4, 0);
        v1[0] = 120;
        v1[3] = 111;
        let mut v2: Vector<SparseStorage<u32>> = Vector::new(4, 0);
        v2[0] = 1;
        v2[2] = 111;
        v2[3] = 1;
        println!("{:?}", v1);
        println!("{:?}", v2);
        // v1 + v2
        let added = &v1 + &v2;
        let muled = &v1 * &v2;
        println!("{:?}", added);
        assert_eq!(added[0], 120 + 1);
        assert_eq!(added[1], 0 + 0);
        assert_eq!(added[2], 0 + 111);
        assert_eq!(added[3], 111 + 1);
        println!("{:?}", muled);
        assert_eq!(muled[0], 120 * 1);
        assert_eq!(muled[1], 0 * 0);
        assert_eq!(muled[2], 0 * 111);
        assert_eq!(muled[3], 111 * 1);
    }
}
