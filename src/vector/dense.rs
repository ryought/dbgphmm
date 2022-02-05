//!
//! Dense storage that uses `std::Vec`
//!
use super::Storage;

/// Dense storage powered by `std::Vec`
///
/// In `DenseStorage`, internal id equals to index.
#[derive(Debug, Clone)]
pub struct DenseStorage<T>(Vec<T>);

impl<T> Storage for DenseStorage<T>
where
    T: Copy,
{
    type Item = T;
    fn new(size: usize, default_value: T) -> DenseStorage<T> {
        DenseStorage(vec![default_value; size])
    }
    #[inline]
    fn size(&self) -> usize {
        self.0.len()
    }
    #[inline]
    fn n_ids(&self) -> usize {
        self.0.len()
    }
    #[inline]
    fn get_by_id(&self, id: usize) -> (usize, T) {
        let index = id;
        let value = self.0[id];
        (index, value)
    }
    #[inline]
    fn get(&self, index: usize) -> &T {
        &self.0[index]
    }
    #[inline]
    fn get_mut(&mut self, index: usize) -> &mut T {
        &mut self.0[index]
    }
}

#[cfg(test)]
mod tests {
    use super::super::Vector;
    use super::*;

    #[test]
    fn dense_storage() {
        // u32
        let mut s: DenseStorage<u32> = DenseStorage::new(10, 0);
        *s.get_mut(5) = 111;
        *s.get_mut(3) = 22;
        assert_eq!(*s.get(0), 0);
        assert_eq!(*s.get(3), 22);
        assert_eq!(*s.get(5), 111);
        assert_eq!(s.size(), 10);
        assert_eq!(s.n_ids(), 10);
        assert_eq!(s.get_by_id(3), (3, 22));

        // clone
        let mut s2 = s.clone();
        assert_eq!(*s2.get(3), 22);
        *s2.get_mut(3) = 21;
        assert_eq!(*s2.get(3), 21);
        assert_eq!(*s.get(3), 22);

        // f64
        let mut s: DenseStorage<f64> = DenseStorage::new(10, 0.0);
        *s.get_mut(5) = 12.11;
        *s.get_mut(3) = 10.0;
        assert_eq!(*s.get(0), 0.0);
        assert_eq!(*s.get(3), 10.0);
        assert_eq!(*s.get(5), 12.11);
    }
    #[test]
    #[should_panic]
    fn dense_storage_outside() {
        let mut s: DenseStorage<u32> = DenseStorage::new(3, 0);
        *s.get_mut(3) = 22;
    }
    #[test]
    fn dense_storage_iter() {
        let mut s: DenseStorage<u32> = DenseStorage::new(4, 5);
        *s.get_mut(0) = 111;
        *s.get_mut(2) = 10;
        let v: Vec<(usize, u32)> = s.iter().collect();
        assert_eq!(v, vec![(0, 111), (1, 5), (2, 10), (3, 5)]);
    }
}
