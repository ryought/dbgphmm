use crate::common::{Sequence, VALID_BASES};
use itertools::Itertools;
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;

///
/// generate random bases of given length from seed
///
pub fn generate(length: usize, seed: u64) -> Sequence {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut seq: Sequence = Vec::with_capacity(length);
    for _ in 0..length {
        // let n: u8 = rng.gen_range(0..4);
        let base = VALID_BASES.choose(&mut rng).unwrap();
        seq.push(*base);
    }
    seq
}

///
/// generate tandem repeat sequence with a unit repeating n_repeat times.
///
pub fn tandem_repeat(unit: &Sequence, n_repeat: usize) -> Sequence {
    let mut seq: Sequence = Vec::with_capacity(unit.len() * n_repeat);
    for _ in 0..n_repeat {
        seq.extend(unit);
    }
    seq
}

/// a edit operation to seq x = x[0], ..., x[n-1].
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Ord, Eq)]
pub struct EditOperation(usize, EditType);

impl EditOperation {
    pub fn new(index: usize, edit_type: EditType) -> EditOperation {
        EditOperation(index, edit_type)
    }
    pub fn index(&self) -> usize {
        self.0
    }
    pub fn edit_type(&self) -> EditType {
        self.1
    }
}

/// a edit operation to seq x = x[0], ..., x[n-1].
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Ord, Eq)]
pub enum EditType {
    /// mutation of x[i] to the specified base
    Mut(u8),
    /// insetion of a base in the left of x[i]
    Ins(u8),
    /// deletion of x[i]
    Del,
}

fn pick_operation<R: Rng>(rng: &mut R, length: usize) -> EditOperation {
    // (1) pick operation
    // (2) pick position
    unimplemented!();
}

fn apply_operations(seq: &mut Sequence, ops: &[EditOperation]) {
    let ops: Vec<_> = ops.iter().sorted().dedup().collect();
    let mut offset: i64 = 0;
    for op in ops {
        let index = op.index() as i64 + offset;
        if index >= 0 && index < seq.len() as i64 {
            let index = index as usize;
            match op.edit_type() {
                EditType::Mut(base) => {
                    seq[index] = base;
                }
                EditType::Ins(base) => {
                    seq.insert(index, base);
                    offset += 1;
                }
                EditType::Del => {
                    seq.remove(index);
                    offset -= 1;
                }
            }
        }
    }
}

///
/// Mutate randomly with
///
/// * Mut(pos)
/// * Ins(pos)
/// * Del(pos)
///
pub fn random_mutation() -> Sequence {
    unimplemented!();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::sequence_to_string;

    #[test]
    fn random_seq_generate() {
        let s = generate(10, 0);
        println!("{:?}", sequence_to_string(&s));
        assert_eq!(s.len(), 10);
        assert_eq!(s, b"CCAATTCACA");

        let s = generate(50, 0);
        println!("{:?}", sequence_to_string(&s));
        assert_eq!(s.len(), 50);
        assert_eq!(s, b"CCAATTCACAAAAACCACACCTTGGCCAAGGTATCGTATCTTGTTGTTGT");

        let s = generate(50, 11);
        println!("{:?}", sequence_to_string(&s));
        assert_eq!(s.len(), 50);
        assert_eq!(s, b"TTGCATCCCATAATTCAGTATAGCCATGGGCTGGGCGCGTAAGGATGCTC");
    }
    #[test]
    fn random_seq_tandem_repeat() {
        let unit = b"ATTT".to_vec();
        let s = tandem_repeat(&unit, 4);
        println!("{}", sequence_to_string(&s));
        assert_eq!(s, b"ATTTATTTATTTATTT");
    }
    #[test]
    fn random_seq_edit_ops() {
        let mut s: Sequence = b"ATCGGATCA".to_vec();
        apply_operations(
            &mut s,
            &[
                EditOperation::new(0, EditType::Del),
                EditOperation::new(0, EditType::Mut(b'X')),
                EditOperation::new(0, EditType::Ins(b'Y')),
                EditOperation::new(1, EditType::Mut(b'X')),
                EditOperation::new(4, EditType::Del),
                EditOperation::new(4, EditType::Mut(b'X')),
                EditOperation::new(4, EditType::Ins(b'Y')),
                EditOperation::new(8, EditType::Del),
                EditOperation::new(8, EditType::Mut(b'X')),
                EditOperation::new(8, EditType::Ins(b'Y')),
            ],
        );
        println!("{}", sequence_to_string(&s));
        assert_eq!(s, b"YXCGYATCY");

        let mut s: Sequence = b"ATCGGATCA".to_vec();
        apply_operations(
            &mut s,
            &[
                EditOperation::new(0, EditType::Mut(b'X')),
                EditOperation::new(2, EditType::Mut(b'X')),
                EditOperation::new(4, EditType::Mut(b'X')),
                EditOperation::new(5, EditType::Mut(b'X')),
                EditOperation::new(8, EditType::Mut(b'X')),
                EditOperation::new(9, EditType::Mut(b'X')),
            ],
        );
        println!("{}", sequence_to_string(&s));
        assert_eq!(s, b"XTXGXXTCX");
    }
}
