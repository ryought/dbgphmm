//!
//! QuadArray
//! store 4 values (0,1,2,3) in u64
//!

/// u64 as [u2; 32]
/// u2 is a pseudo-type with 4-values
#[derive(PartialEq, PartialOrd, Eq, Hash, Clone, Copy)]
pub struct QuadArray(pub u64);

impl QuadArray {
    pub fn new() -> QuadArray {
        QuadArray(0)
    }
    pub fn empty() -> QuadArray {
        QuadArray(0)
    }
    /// get QuadArray[index]
    pub fn get(&self, index: usize) -> u64 {
        assert!(index >= 0);
        assert!(index < 32);
        (self.0 >> (index * 2)) & 0b11u64
    }
    /// set QuadArray[index]
    pub fn set(&mut self, index: usize, value: u64) {
        assert!(index >= 0);
        assert!(index < 32);
        assert!(value >= 0);
        assert!(value < 4);
        let shift = index * 2;
        self.0 = self.0 & !(0b11u64 << shift) | ((value as u64) << shift)
    }
    /// shift all elements to front
    /// [0,3,2,3,1] -> [3,2,3,1,0]
    pub fn shift_front(&mut self) {
        self.0 = self.0 >> 2
    }
    /// shift all elements to end
    /// [0,3,2,3,1] -> [0,0,3,2,3]
    pub fn shift_back(&mut self) {
        self.0 = self.0 << 2
    }
    /// Convert to [u8; 32]
    fn to_array(&self) -> [u8; 32] {
        let mut arr = [0; 32];
        for i in 0..32 {
            arr[i] = self.get(i) as u8;
        }
        arr
    }
}

//
// Tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quadarray() {
        let mut q = QuadArray(0b11001011_u64);
        assert_eq!(q.get(0), 0b11);
        assert_eq!(q.get(1), 0b10);
        assert_eq!(q.get(2), 0b00);
        assert_eq!(q.get(3), 0b11);
        assert_eq!(q.get(4), 0b00);
        q.set(1, 0b11);
        assert_eq!(q.get(0), 0b11);
        assert_eq!(q.get(1), 0b11);
        assert_eq!(q.get(2), 0b00);
        assert_eq!(q.get(3), 0b11);
        assert_eq!(q.get(4), 0b00);
        assert_eq!(
            q.to_array(),
            [
                3, 3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ]
        );
        q.shift_front();
        assert_eq!(q.get(0), 0b11);
        assert_eq!(q.get(1), 0b00);
        assert_eq!(q.get(2), 0b11);
        assert_eq!(q.get(3), 0b00);
        assert_eq!(q.get(4), 0b00);
    }
}
