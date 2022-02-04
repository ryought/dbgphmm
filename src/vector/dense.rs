//!
//! Dense storage that uses `std::Vec`
//!
use super::{IterableStorage, Storage, Vector};

/// Dense storage powered by `std::Vec`
#[derive(Debug, Clone)]
struct DenseStorage<T>(Vec<T>);

impl<T> Storage for DenseStorage<T>
where
    T: Copy,
{
    type Item = T;
    fn new(size: usize, default_value: T) -> DenseStorage<T> {
        DenseStorage(vec![default_value; size])
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

impl<'a, T> IterableStorage<'a> for DenseStorage<T>
where
    T: Copy + 'a,
{
    type IndexIterator = DenseStorageIterator<'a, T>;
    fn indexiter(&'a self) -> Self::IndexIterator {
        DenseStorageIterator {
            index: 0,
            storage: &self.0,
        }
    }
}

/// Iterator on DenseStorage (index, item)
struct DenseStorageIterator<'a, T: Copy + 'a> {
    /// current index on the storage
    index: usize,
    /// reference to the original storage
    storage: &'a [T],
}

impl<'a, T: Copy + 'a> Iterator for DenseStorageIterator<'a, T> {
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
        let mut s: DenseStorage<u32> = DenseStorage::new(10, 0);
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
        let v: Vec<(usize, u32)> = s.indexiter().collect();
        assert_eq!(v, vec![(0, 111), (1, 5), (2, 10), (3, 5)]);
    }
    #[test]
    fn dense_storage_vector() {
        let mut v: Vector<DenseStorage<u32>> = Vector::new(10, 0);
        v[0] = 100;
        v[3] = 222;
        assert_eq!(v[0], 100);
        assert_eq!(v[1], 0);
        let w: Vec<(usize, u32)> = v.iter().collect();
        println!("{:?}", v);
        println!("{:?}", w);
        assert_eq!(
            w,
            vec![
                (0, 100),
                (1, 0),
                (2, 0),
                (3, 222),
                (4, 0),
                (5, 0),
                (6, 0),
                (7, 0),
                (8, 0),
                (9, 0)
            ]
        );
    }
    #[test]
    fn dense_storage_vector_add() {
        let mut v1: Vector<DenseStorage<u32>> = Vector::new(4, 0);
        v1[0] = 120;
        v1[3] = 111;
        let mut v2: Vector<DenseStorage<u32>> = Vector::new(4, 0);
        v2[0] = 1;
        v2[2] = 111;
        v2[3] = 1;
        println!("{:?}", v1);
        println!("{:?}", v2);
        // v1 + v2
        println!("{:?}", &v1 + &v2);
    }
}
