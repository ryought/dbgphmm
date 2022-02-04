//!
//! Dense storage that uses `std::Vec`
//!
use super::{IterableStorage, Storage};

/// Dense storage
#[derive(Debug, Clone)]
struct VectorStorage<T>(Vec<T>);

impl<T> Storage for VectorStorage<T>
where
    T: Copy,
{
    type Item = T;
    fn new(size: usize, default_value: T) -> VectorStorage<T> {
        VectorStorage(vec![default_value; size])
    }
    fn size(&self) -> usize {
        self.0.len()
    }
    fn get(&self, index: usize) -> &T {
        &self.0[index]
    }
    fn get_mut(&mut self, index: usize) -> &mut T {
        &mut self.0[index]
    }
}

impl<'a, T> IterableStorage<'a> for VectorStorage<T>
where
    T: Copy + 'a,
{
    type IndexIterator = VectorStorageIterator<'a, T>;
    fn indexiter(&'a self) -> Self::IndexIterator {
        VectorStorageIterator {
            index: 0,
            storage: &self.0,
        }
    }
}

/// Iterator on VectorStorage
struct VectorStorageIterator<'a, T: Copy + 'a> {
    /// current index on the vector
    index: usize,
    /// reference to the original storage
    storage: &'a [T],
}

impl<'a, T: Copy + 'a> Iterator for VectorStorageIterator<'a, T> {
    type Item = (usize, T);
    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.storage.len() {
            let index = self.index;
            let value = self.storage[index];
            self.index += 1;
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
    fn dense_storage() {
        // u32
        let mut s: VectorStorage<u32> = VectorStorage::new(10, 0);
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
        let mut s: VectorStorage<f64> = VectorStorage::new(10, 0.0);
        *s.get_mut(5) = 12.11;
        *s.get_mut(3) = 10.0;
        assert_eq!(*s.get(0), 0.0);
        assert_eq!(*s.get(3), 10.0);
        assert_eq!(*s.get(5), 12.11);
    }
    #[test]
    #[should_panic]
    fn dense_storage_outside() {
        let mut s: VectorStorage<u32> = VectorStorage::new(3, 0);
        *s.get_mut(3) = 22;
    }
    #[test]
    fn dense_storage_iter() {
        let mut s: VectorStorage<u32> = VectorStorage::new(4, 5);
        *s.get_mut(0) = 111;
        *s.get_mut(2) = 10;
        let v: Vec<(usize, u32)> = s.indexiter().collect();
        assert_eq!(v, vec![(0, 111), (1, 5), (2, 10), (3, 5)]);
    }
}
