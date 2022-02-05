//!
//! `Vector` Wrapper of fixed size table
//!
//!
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign};
pub mod dense;
pub mod graph;
pub mod index;
// pub mod sparse;
// pub use sparse::SparseStorage;
pub use dense::DenseStorage;
pub use graph::NodeVec;
pub use petgraph::graph::IndexType;
use std::marker::PhantomData;

/// Backend storage of `Vector`
/// an abstruction of a vec with fixed size that is readable/writable
/// by index.
///
/// * `new`
///     create storage with fixed size and filled with the default value
/// * `size`
///     get the fixed size
/// * `get`
///     get the reference to the value in the index
/// * `get_mut`
///     get the mutable reference to the value in the index
/// * `iter`
///     get an iterator of (index, value) whose value is not default.
///
/// Internaly, mapping of index (0 <= index < size) to value
/// is stored with internal element id. (0 <= id < # of stored elements)
/// (index) -> (id) -> (value)
///
/// * `get_id`
///     conversion from index to id
/// * `get_by_id`
///     get reference by id
/// * `get_mut_by_id`
///     get mutable reference by id
///
pub trait Storage: Clone + Sized {
    /// Item type that this storage stores.
    ///
    type Item: Copy;
    ///
    /// Create a new storage with fixed size and filled with the default value
    fn new(size: usize, default_value: Self::Item) -> Self;
    ///
    /// Get the size of this storage
    fn size(&self) -> usize;
    ///
    /// Get internal id of the index
    fn get_id(&self, index: usize) -> usize;
    ///
    /// Get the reference to the value at the internal id
    fn get_by_id(&self, id: usize) -> &Self::Item;
    ///
    /// Get the mutable reference at the internal id
    fn get_mut_by_id(&mut self, id: usize) -> &mut Self::Item;
    ///
    /// Get the reference to the value at the given index
    fn get(&self, index: usize) -> &Self::Item {
        let id = self.get_id(index);
        self.get_by_id(id)
    }
    ///
    /// Get the mutable reference to the given index
    fn get_mut(&mut self, index: usize) -> &mut Self::Item {
        let id = self.get_id(index);
        self.get_mut_by_id(id)
    }
    ///
    /// get an iterator of (usize, Self::Item) on the storage
    fn iter<'a>(&'a self) -> StorageIterator<'a, Self> {
        StorageIterator {
            id: 0,
            storage: self,
        }
    }
}

///
/// Iterator struct of Storage
pub struct StorageIterator<'a, S: Storage> {
    /// current internal id
    id: usize,
    /// reference of the storage
    storage: &'a S,
}

impl<'a, S: Storage> Iterator for StorageIterator<'a, S> {
    type Item = (usize, S::Item);
    fn next(&mut self) -> Option<Self::Item> {
        let ret = if self.id == 0 {
            // TODO
            Some((0, *self.storage.get(0)))
        } else {
            None
        };
        self.id += 1;
        ret
    }
}

/// `Vector` struct
/// It generalized of (1) its item and (2) its backend storage.
#[derive(Clone, Debug)]
pub struct Vector<S: Storage, Ix: IndexType = usize> {
    /// Backend storage of the Vector
    storage: S,
    ty: PhantomData<Ix>,
}

impl<S: Storage> Vector<S> {
    /// Create a new Vector, with fixed size and filled by default_value.
    pub fn new(size: usize, default_value: S::Item) -> Vector<S> {
        Vector {
            storage: S::new(size, default_value),
            ty: PhantomData,
        }
    }
    /// Get an (virtual) size of the storage
    pub fn len(&self) -> usize {
        self.storage.size()
    }
    /// Get an iterator on (index, item).
    pub fn iter<'a>(&'a self) -> StorageIterator<'a, S> {
        self.storage.iter()
    }
}

/// Implement index access, vec[i]
impl<S: Storage> Index<usize> for Vector<S> {
    type Output = S::Item;
    fn index(&self, index: usize) -> &Self::Output {
        self.storage.get(index)
    }
}

/// Implement index write access, vec[i] = 10
impl<S: Storage> IndexMut<usize> for Vector<S> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.storage.get_mut(index)
    }
}

/*
/// Implement addition `+` between two vecs
/// if the item of vec supports addition
impl<'a, 'b, S> Add<&'a Vector<S>> for &'b Vector<S>
where
    S: Storage,
    S::Item: Add<Output = S::Item>,
{
    type Output = Vector<S>;
    fn add(self, other: &'a Vector<S>) -> Self::Output {
        assert_eq!(self.len(), other.len());
        let mut ret = self.clone();
        for (index, value) in other.iter() {
            let old_value = ret[index];
            ret[index] = old_value + value;
        }
        ret
    }
}

/// Implement addition with assignment `+=` between two vecs
/// if the item of vec supports addition
/// This does not cause re-allocation
impl<'a, S> AddAssign<&'a Vector<S>> for Vector<S>
where
    S: Storage,
    S::Item: Add<Output = S::Item>,
{
    fn add_assign(&mut self, other: &'a Vector<S>) {
        assert_eq!(self.len(), other.len());
        for (index, value) in other.iter() {
            self[index] = self[index] + value;
        }
    }
}

/// Implement multiplication `*` between two vecs
/// if the item of vec supports multiplication
impl<'a, 'b, S> Mul<&'a Vector<S>> for &'b Vector<S>
where
    S: IterableStorage<'a>,
    S::Item: Mul<Output = S::Item>,
{
    type Output = Vector<S>;
    fn mul(self, other: &'a Vector<S>) -> Self::Output {
        assert_eq!(self.len(), other.len());
        let mut ret = self.clone();
        for (index, value) in other.iter() {
            let old_value = ret[index];
            ret[index] = old_value * value;
        }
        ret
    }
}

/// Implement multiplication with assignment `*=` between two vecs
/// if the item of vec supports multiplication
/// This does not cause re-allocation
impl<'a, S> MulAssign<&'a Vector<S>> for Vector<S>
where
    S: IterableStorage<'a>,
    S::Item: Mul<Output = S::Item>,
{
    fn mul_assign(&mut self, other: &'a Vector<S>) {
        assert_eq!(self.len(), other.len());
        for (index, value) in other.iter() {
            self[index] = self[index] * value;
        }
    }
}
*/

#[cfg(test)]
mod tests {
    use super::*;
}
