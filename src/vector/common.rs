//!
//! `vector` for
//!

use std::ops::{Add, Index, IndexMut};

///
trait Storage: Clone {
    type Item: Copy;
    fn new(size: usize, default_value: Self::Item) -> Self;
    fn size(&self) -> usize;
    fn get(&self, index: usize) -> &Self::Item;
    fn get_mut(&mut self, index: usize) -> &mut Self::Item;
}

trait IterableStorage<'a>: Storage {
    type IndexIterator: Iterator<Item = (usize, Self::Item)>;
    fn indexiter(&'a self) -> Self::IndexIterator;
}
