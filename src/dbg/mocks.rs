use super::hashdbg_v2::HashDbg;
use super::impls::SimpleDbg;
use crate::kmer::veckmer::VecKmer;

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
