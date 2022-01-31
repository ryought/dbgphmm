use dbgphmm::kmer::common::{KmerBase, KmerLike, NullableKmer};
use dbgphmm::kmer::tinykmer::TinyKmer;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let k: usize = args[1].to_string().parse().unwrap();
    let kmer = TinyKmer::from(b"ATCGGGAT");
    println!("k={} kmer={}", k, kmer);
}
