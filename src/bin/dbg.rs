use dbgphmm::dbg::hashdbg_v2::HashDbg;
use dbgphmm::dbg::impls::SimpleDbg;
use dbgphmm::kmer::veckmer::VecKmer;

fn main() {
    // let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"AAAGCTTGATT");
    let mut hd: HashDbg<VecKmer> = HashDbg::new(4);
    hd.add_seq(b"AAAGCTTGATT");
    hd.add_seq(b"CGTATC");
    let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_hashdbg(&hd);
    println!("{}", dbg);
}
