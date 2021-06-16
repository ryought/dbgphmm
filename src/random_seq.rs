use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;

pub fn generate(length: usize, seed: u64) -> Vec<u8> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut seq: Vec<u8> = Vec::with_capacity(length);
    for _ in 0..length {
        // let n: u8 = rng.gen_range(0..4);
        let bases = [b'A', b'C', b'G', b'T'];
        let base = bases.choose(&mut rng).unwrap();
        seq.push(*base);
    }
    seq
}
