use crate::dbg::{DbgHash, DBG};
use crate::kmer::kmer::Kmer;
use bio::io::fasta;
use log::warn;

pub fn sanitize_bases(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .enumerate()
        .map(|(i, base)| match base {
            b'A' | b'a' => b'A',
            b'C' | b'c' => b'C',
            b'G' | b'g' => b'G',
            b'T' | b't' => b'T',
            b'N' | b'n' => {
                warn!("ambiguous detected `n` in bases[{}]", i);
                b'N'
            }
            &c => {
                warn!("informal base `{}` detected in bases[{}]", c as char, i);
                b'N'
            }
        })
        .collect()
}

pub fn parse_kmers_and_copy_nums(filename: &str, k: usize) -> (Vec<Kmer>, Vec<u32>) {
    // and returns (kmers, copy_nums) that can be used as an input
    // of DbgPHMM::new.
    let reader = fasta::Reader::from_file(filename).unwrap();
    let mut d = DbgHash::new();
    for result in reader.records() {
        let record = result.unwrap();
        for window in sanitize_bases(record.seq()).windows(k) {
            let kmer = Kmer::from(window);
            d.add(kmer, 1);
        }
    }
    let kmers = d.kmers();
    let copy_nums: Vec<u32> = kmers.iter().map(|kmer| d.find(kmer)).collect();
    // println!("{}", d.as_dot());
    (kmers, copy_nums)
}

pub fn parse_reads(filename: &str) -> Vec<Vec<u8>> {
    let mut reader = fasta::Reader::from_file(filename).unwrap();
    let mut reads: Vec<Vec<u8>> = Vec::new();
    for result in reader.records() {
        let record = result.unwrap();
        reads.push(sanitize_bases(record.seq()));
    }
    reads
}

pub fn read2() {
    let s = b"ATCGATTCGATCGATTCGAT";
    println!("all {:?}", s);
    for w in s.windows(5) {
        println!("window {:?}", w);
    }
}
