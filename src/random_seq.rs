use crate::common::{different_bases, Sequence, VALID_BASES};
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
        let base = pick_random_base(&mut rng);
        seq.push(base);
    }
    seq
}

///
/// concatenete two sequences
///
pub fn join(a: Sequence, b: Sequence) -> Sequence {
    [a, b].concat()
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

fn pick_different_base<R: Rng>(rng: &mut R, base: u8) -> u8 {
    *different_bases(base).choose(rng).unwrap()
}

fn pick_random_base<R: Rng>(rng: &mut R) -> u8 {
    *VALID_BASES.choose(rng).unwrap()
}

fn pick_operation<R: Rng>(rng: &mut R, seq: &Sequence, profile: MutationProfile) -> EditOperation {
    // (1) pick position
    let index: usize = rng.gen_range(0..seq.len());
    let base = seq[index];
    // (2) pick edit type
    let base_mut = pick_different_base(rng, base);
    let base_ins = pick_random_base(rng);
    let edit_type = *[
        EditType::Mut(base_mut),
        EditType::Ins(base_ins),
        EditType::Del,
    ]
    .choose_weighted(rng, |edit_type| match edit_type {
        EditType::Mut(_) => profile.p_mut,
        EditType::Ins(_) => profile.p_ins,
        EditType::Del => profile.p_del,
    })
    .unwrap();

    EditOperation::new(index, edit_type)
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

#[derive(Clone, Debug, Copy)]
pub struct MutationProfile {
    /// divergence rate
    divergence_rate: f64,
    /// probability of mutation
    p_mut: f64,
    /// probability of insertion
    p_ins: f64,
    /// probability of deletion
    p_del: f64,
}

impl MutationProfile {
    pub fn uniform(divergence_rate: f64) -> MutationProfile {
        MutationProfile {
            divergence_rate,
            p_mut: 1.0,
            p_ins: 1.0,
            p_del: 1.0,
        }
    }
}

///
/// Mutate randomly with MutationProfile
///
pub fn random_mutation(
    seq: &Sequence,
    profile: MutationProfile,
    seed: u64,
) -> (Sequence, Vec<EditOperation>) {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let n_mutations = (seq.len() as f64 * profile.divergence_rate)
        .round()
        .max(0.0) as usize;

    // create operations
    let ops: Vec<_> = (0..n_mutations)
        .map(|_| pick_operation(&mut rng, seq, profile))
        .collect();

    // apply the operations to copied seq
    let mut ret = seq.to_vec();
    apply_operations(&mut ret, &ops);
    (ret, ops)
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
    fn random_seq_join() {
        let a = b"ATTT".to_vec();
        let b = b"GGGG".to_vec();
        let c = join(a, b);
        println!("{}", sequence_to_string(&c));
        assert_eq!(c, b"ATTTGGGG");
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
    #[test]
    fn random_seq_pick_operation() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
        let s = b"ATCGATCGT".to_vec();
        let ops: Vec<_> = (0..10)
            .map(|_| {
                pick_operation(
                    &mut rng,
                    &s,
                    MutationProfile {
                        divergence_rate: 0.01,
                        p_mut: 1.0,
                        p_ins: 1.0,
                        p_del: 1.0,
                    },
                )
            })
            .collect();
        println!("{:?}", ops);
        assert_eq!(
            ops,
            vec![
                EditOperation::new(3, EditType::Ins(65)),
                EditOperation::new(0, EditType::Mut(84)),
                EditOperation::new(1, EditType::Mut(67)),
                EditOperation::new(0, EditType::Mut(71)),
                EditOperation::new(3, EditType::Del),
                EditOperation::new(1, EditType::Mut(71)),
                EditOperation::new(3, EditType::Mut(65)),
                EditOperation::new(1, EditType::Ins(71)),
                EditOperation::new(7, EditType::Del),
                EditOperation::new(4, EditType::Ins(65))
            ]
        );
    }
    #[test]
    fn random_seq_random_mutation() {
        let s = generate(50, 0);

        // 1
        let (t, ops) = random_mutation(
            &s,
            MutationProfile {
                divergence_rate: 0.0001,
                p_mut: 1.0,
                p_ins: 1.0,
                p_del: 1.0,
            },
            0,
        );
        println!("{}", sequence_to_string(&t));
        println!("{:?}", ops);
        assert_eq!(ops.len(), 0);
        assert_eq!(s, t);
        assert_eq!(t, b"CCAATTCACAAAAACCACACCTTGGCCAAGGTATCGTATCTTGTTGTTGT");

        // 2
        let (t, ops) = random_mutation(
            &s,
            MutationProfile {
                divergence_rate: 0.1,
                p_mut: 1.0,
                p_ins: 0.0,
                p_del: 0.0,
            },
            0,
        );
        println!("{}", sequence_to_string(&s));
        println!("{}", sequence_to_string(&t));
        println!("{:?}", ops);
        // all ops must be Mut
        assert!(ops.iter().all(|op| {
            if let EditType::Mut(_) = op.edit_type() {
                true
            } else {
                false
            }
        }));
        assert_eq!(t, b"CCAATACACAAAAAACGCACCTTGACCAAGGAATCGTATCTTGTTGTTGT");

        // 3
        let (t, ops) = random_mutation(
            &s,
            MutationProfile {
                divergence_rate: 0.2,
                p_mut: 1.0,
                p_ins: 1.0,
                p_del: 1.0,
            },
            0,
        );
        println!("{}", sequence_to_string(&s));
        println!("{}", sequence_to_string(&t));
        assert_eq!(t, b"CCAATATCCACAGAAAACGAACTTGCCAAGGCTTTCGTATCTTGTTGTTGT");
    }
}
