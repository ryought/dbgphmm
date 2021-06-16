use crate::dbg::{DbgHash, DBG};
use crate::kmer::kmer::Kmer;
use bio::io::fasta;

pub fn parse_kmers_and_copy_nums(filename: &str, k: usize) -> (Vec<Kmer>, Vec<u32>) {
    // and returns (kmers, copy_nums) that can be used as an input
    // of DbgPHMM::new.
    let mut reader = fasta::Reader::from_file(filename).unwrap();
    let mut d = DbgHash::new();
    for result in reader.records() {
        let record = result.unwrap();
        for window in record.seq().windows(k) {
            let kmer = Kmer::from(window);
            d.add(kmer, 1);
        }
    }
    let kmers = d.kmers();
    let copy_nums: Vec<u32> = kmers.iter().map(|kmer| d.find(kmer)).collect();
    eprintln!("is_consistent: {}", d.is_copy_number_consistent());
    println!("{}", d.as_dot());
    (kmers, copy_nums)
}

pub fn parse_reads(filename: &str) -> Vec<Vec<u8>> {
    let mut reader = fasta::Reader::from_file(filename).unwrap();
    let mut reads: Vec<Vec<u8>> = Vec::new();
    for result in reader.records() {
        let record = result.unwrap();
        reads.push(record.seq().to_vec());
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
