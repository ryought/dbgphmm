use arrayvec::ArrayVec;
use std::ops::{Index, IndexMut};

#[derive(Debug)]
pub struct SparseVec<T: Copy> {
    size: usize,
    index: ArrayVec<usize, 10>,
    value: ArrayVec<T, 10>,
}

impl<T: Copy> SparseVec<T> {
    pub fn new(size: usize) -> SparseVec<T> {
        SparseVec {
            size,
            index: ArrayVec::<usize, 10>::new(),
            value: ArrayVec::<T, 10>::new(),
        }
    }
    pub fn read(&self, i: usize) -> Option<T> {
        for j in 0..self.index.len() {
            if self.index[j] == i {
                match self.value.get(j) {
                    Some(x) => return Some(x.clone()),
                    None => return None,
                }
            }
        }
        return None;
    }
    pub fn write(&mut self, i: usize, x: T) {
        self.index.push(i);
        self.value.push(x);
    }
}

impl<T> Index<usize> for SparseVec<T> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        println!("acccessing {:?}", index);
        match index {
            Side::Left => &self.left,
            Side::Right => &self.right,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_sparse() {
        let mut v = SparseVec::<usize>::new(10);
        v.write(10, 1918);
        v.write(0, 18);
        println!("{:?}", v);
        println!("{:?}", v.read(10));
        println!("{:?}", v.read(20));
        println!("{:?}", v.read(0));
    }
}
