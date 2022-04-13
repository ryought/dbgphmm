use super::hashdbg_v2::HashDbg;
use super::impls::SimpleDbg;
use crate::common::Sequence;
use crate::kmer::veckmer::VecKmer;
use crate::random_seq::generate;

pub fn mock_base() -> SimpleDbg<VecKmer> {
    let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"ATCGGCT");
    SimpleDbg::from_hashdbg(&hd)
}

pub fn mock_simple() -> SimpleDbg<VecKmer> {
    let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"AAAGCTTGATT");
    SimpleDbg::from_hashdbg(&hd)
}

pub fn mock_two_seqs() -> SimpleDbg<VecKmer> {
    let mut hd: HashDbg<VecKmer> = HashDbg::new(4);
    hd.add_seq(b"AAAGCTTGATT");
    hd.add_seq(b"CGTATC");
    SimpleDbg::from_hashdbg(&hd)
}

pub fn mock_rep() -> SimpleDbg<VecKmer> {
    let mut hd: HashDbg<VecKmer> = HashDbg::new(4);
    hd.add_seq(b"AAAAAAAAAAAAA");
    hd.add_seq(b"CCCCCCCCCCCCCC");
    SimpleDbg::from_hashdbg(&hd)
}

///
/// AACTAGCTT x1
/// CCGTAGGGC x1
///
/// `TAG` is intersecting k-1mer.
///
pub fn mock_intersection() -> SimpleDbg<VecKmer> {
    let mut hd: HashDbg<VecKmer> = HashDbg::new(4);
    hd.add_seq(b"AACTAGCTT");
    hd.add_seq(b"CCGTAGGGC");
    SimpleDbg::from_hashdbg(&hd)
}

///
/// ATAGCT x1
/// TTAGAT x1
///
/// `TAG` is intersecting k-1-mer.
///
pub fn mock_intersection_small() -> SimpleDbg<VecKmer> {
    let mut hd: HashDbg<VecKmer> = HashDbg::new(4);
    hd.add_seq(b"ATAGCT");
    hd.add_seq(b"TAGGAT");
    SimpleDbg::from_hashdbg(&hd)
}

///
/// crate a mock `SimpleDbg` from single random sequence of given length
///
pub fn mock_random_with_seq(k: usize, length: usize) -> (SimpleDbg<VecKmer>, Sequence) {
    let seq = generate(length, 1);
    let mut hd: HashDbg<VecKmer> = HashDbg::new(k);
    hd.add_seq(&seq);
    (SimpleDbg::from_hashdbg(&hd), seq)
}

pub fn mock_random(k: usize, length: usize) -> SimpleDbg<VecKmer> {
    let (dbg, _) = mock_random_with_seq(k, length);
    dbg
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dbg_mock_random() {
        let dbg = mock_random(8, 100);
        println!("{}", dbg);
    }
}
