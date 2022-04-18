use crate::common::Sequence;
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
        let bases = [b'A', b'C', b'G', b'T'];
        let base = bases.choose(&mut rng).unwrap();
        seq.push(*base);
    }
    seq
}

pub fn random_mutation() -> Sequence {
    unimplemented!();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::sequence_to_string;

    #[test]
    fn random_seq() {
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
}
