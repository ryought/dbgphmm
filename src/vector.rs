//!
//! `Vector` Wrapper of fixed size table
//!
//!
use std::ops::{Add, AddAssign, Div, Index, IndexMut, Mul, MulAssign};
pub mod dense;
pub mod graph;
pub mod index;
pub mod sparse;
pub mod test;
pub use dense::DenseStorage;
pub use graph::{EdgeVec, NodeVec};
pub use index::Indexable;
pub use sparse::SparseStorage;
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
/// Internaly, mapping of index (0 <= index < self.size()) to value
/// is stored with internal element id. (0 <= id < self.n_ids())
/// (index) -> (id) -> (value)
///
pub trait Storage: Clone + Sized {
    /// Item type that this storage stores.
    ///
    type Item: Copy + PartialEq;
    ///
    /// Create a new storage with fixed size and filled with the default value
    fn new(size: usize, default_value: Self::Item) -> Self;
    ///
    /// Get the size of this storage
    fn size(&self) -> usize;
    ///
    /// Get the number of elements.
    /// `0 <= id < self.n_ids()` should be satisfied.
    fn n_ids(&self) -> usize;
    ///
    /// Get the reference to the value at the internal id
    fn get_by_id(&self, id: usize) -> (usize, Self::Item);
    ///
    /// Get the reference to the value at the given index
    fn get(&self, index: usize) -> &Self::Item;
    ///
    /// Get the mutable reference to the given index
    fn get_mut(&mut self, index: usize) -> &mut Self::Item;
    ///
    /// get an iterator of (usize, Self::Item) on the storage
    fn iter<'a>(&'a self) -> StorageIterator<'a, Self> {
        StorageIterator {
            id: 0,
            storage: self,
        }
    }
    ///
    /// Convert to the DenseStorage with same contents
    fn to_dense(&self) -> DenseStorage<Self::Item>;
    ///
    /// Convert to the SparseStorage with same contents
    /// with specifying `default_value` in SparseStorage.
    fn to_sparse(&self, default_value: Self::Item) -> SparseStorage<Self::Item>;
    ///
    /// Check if this is dense storage or not
    fn is_dense() -> bool;
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
        if self.id < self.storage.n_ids() {
            let id = self.id;
            let item = self.storage.get_by_id(id);
            self.id += 1;
            Some(item)
        } else {
            None
        }
    }
}

/// `Vector` struct
///
/// It generalized of
///
/// 1. item type `Storage::Item`
/// 2. backend storage `S: Storage`
/// 3. index type `Ix: Indexable`
///
#[derive(Clone, Debug)]
pub struct Vector<S: Storage, Ix: Indexable = usize> {
    /// Backend storage of the Vector
    storage: S,
    /// Hidden marker of index type
    ty: PhantomData<Ix>,
}

impl<S: Storage, Ix: Indexable> Vector<S, Ix> {
    /// Create a new Vector, with fixed size and filled by default_value.
    pub fn new(size: usize, default_value: S::Item) -> Vector<S, Ix> {
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
    pub fn iter<'a>(&'a self) -> impl 'a + Iterator<Item = (Ix, S::Item)> {
        self.storage.iter().map(|(i, v)| (Ix::new(i), v))
    }
    /// Convert to the DenseStorage-backed vector.
    pub fn to_dense(&self) -> Vector<DenseStorage<S::Item>, Ix> {
        Vector {
            storage: self.storage.to_dense(),
            ty: PhantomData,
        }
    }
    /// Convert to the SparseStorage-backed vector.
    pub fn to_sparse(&self, default_value: S::Item) -> Vector<SparseStorage<S::Item>, Ix> {
        Vector {
            storage: self.storage.to_sparse(default_value),
            ty: PhantomData,
        }
    }
    pub fn is_dense(&self) -> bool {
        S::is_dense()
    }
}

/// Implement index access, vec[i]
impl<S: Storage, Ix: Indexable> Index<Ix> for Vector<S, Ix> {
    type Output = S::Item;
    fn index(&self, index: Ix) -> &Self::Output {
        self.storage.get(index.index())
    }
}

/// Implement index write access, vec[i] = 10
impl<S: Storage, Ix: Indexable> IndexMut<Ix> for Vector<S, Ix> {
    fn index_mut(&mut self, index: Ix) -> &mut Self::Output {
        self.storage.get_mut(index.index())
    }
}

/// Implement addition `+` between two vecs
/// if the item of vec supports addition
impl<'a, 'b, Sa, Sb, Ix> Add<&'a Vector<Sa, Ix>> for &'b Vector<Sb, Ix>
where
    Sa: Storage,
    Sb: Storage<Item = Sa::Item>,
    Sa::Item: Add<Output = Sa::Item>,
    Ix: Indexable,
{
    type Output = Vector<Sb, Ix>;
    fn add(self, other: &'a Vector<Sa, Ix>) -> Self::Output {
        assert_eq!(self.len(), other.len());
        let mut ret = self.clone();
        for (index, value) in other.iter() {
            let old_value = ret[index];
            ret[index] = old_value + value;
        }
        ret
    }
}

/// add constant to vector
/// `Vector<S> + S::Item = Vector<S>`
///
/// TODO
/// does not calculate the correct values for `SparseStorage`
/// because it cannot modify the `default_value`.
impl<'a, 'b, S, Ix> Add<S::Item> for Vector<S, Ix>
where
    S: Storage,
    S::Item: Add<Output = S::Item> + Copy,
    Ix: Indexable,
{
    type Output = Vector<S, Ix>;
    fn add(mut self, other: S::Item) -> Self::Output {
        // TODO if sparse, default value should be modified
        // currently, Add<S::Item> can only be used with Vector<Dense>.
        let n = self.storage.n_ids();
        for id in 0..n {
            let (index, value) = self.storage.get_by_id(id);
            *self.storage.get_mut(index) = value + other;
        }
        self
    }
}

/// Implement addition with assignment `+=` between two vecs
/// if the item of vec supports addition
/// This does not cause re-allocation
impl<'a, S, Ix> AddAssign<&'a Vector<S, Ix>> for Vector<S, Ix>
where
    S: Storage,
    S::Item: Add<Output = S::Item>,
    Ix: Indexable,
{
    fn add_assign(&mut self, other: &'a Vector<S, Ix>) {
        assert_eq!(self.len(), other.len());
        for (index, value) in other.iter() {
            self[index] = self[index] + value;
        }
    }
}

/// Implement multiplication `*` between two vecs
/// if the item of vec supports multiplication
impl<'a, 'b, Sa, Sb, Ix> Mul<&'a Vector<Sa, Ix>> for &'b Vector<Sb, Ix>
where
    Sa: Storage,
    Sb: Storage<Item = Sa::Item>,
    Sa::Item: Mul<Output = Sa::Item>,
    Ix: Indexable,
{
    type Output = Vector<Sb, Ix>;
    fn mul(self, other: &'a Vector<Sa, Ix>) -> Self::Output {
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
impl<'a, S, Ix> MulAssign<&'a Vector<S, Ix>> for Vector<S, Ix>
where
    S: Storage,
    S::Item: Mul<Output = S::Item>,
    Ix: Indexable,
{
    fn mul_assign(&mut self, other: &'a Vector<S, Ix>) {
        assert_eq!(self.len(), other.len());
        for (index, value) in other.iter() {
            self[index] = self[index] * value;
        }
    }
}

/// multiply a constant to vector
/// `Vector<S> * S::Item = Vector<S>`
///
/// TODO
/// does not calculate the correct values for `SparseStorage`
/// because it cannot modify the `default_value`.
impl<'a, 'b, S, Ix> Mul<S::Item> for Vector<S, Ix>
where
    S: Storage,
    S::Item: Mul<Output = S::Item> + Copy,
    Ix: Indexable,
{
    type Output = Vector<S, Ix>;
    fn mul(mut self, other: S::Item) -> Self::Output {
        // TODO if sparse, default value should be modified
        // currently, this can only be used with Vector<Dense>.
        let n = self.storage.n_ids();
        for id in 0..n {
            let (index, value) = self.storage.get_by_id(id);
            *self.storage.get_mut(index) = value * other;
        }
        self
    }
}

/// divide-by a constant to vector
/// `Vector<S> / S::Item = Vector<S>`
///
/// TODO
/// does not calculate the correct values for `SparseStorage`
/// because it cannot modify the `default_value`.
impl<'a, 'b, S, Ix> Div<S::Item> for Vector<S, Ix>
where
    S: Storage,
    S::Item: Div<Output = S::Item> + Copy,
    Ix: Indexable,
{
    type Output = Vector<S, Ix>;
    fn div(mut self, other: S::Item) -> Self::Output {
        // TODO if sparse, default value should be modified
        // currently, this can only be used with Vector<Dense>.
        let n = self.storage.n_ids();
        for id in 0..n {
            let (index, value) = self.storage.get_by_id(id);
            *self.storage.get_mut(index) = value / other;
        }
        self
    }
}
