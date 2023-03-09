use crate::common::NULL_BASE;
use crate::kmer::kmer::{Kmer, KmerLike};
use bio::io::fasta;
use log::warn;
use std::io;

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
                NULL_BASE
            }
            &c => {
                warn!("informal base `{}` detected in bases[{}]", c as char, i);
                NULL_BASE
            }
        })
        .collect()
}

pub fn parse_seqs(filename: &str) -> Vec<Vec<u8>> {
    let reader = fasta::Reader::from_file(filename).unwrap();
    let mut reads: Vec<Vec<u8>> = Vec::new();
    for result in reader.records() {
        let record = result.unwrap();
        reads.push(sanitize_bases(record.seq()));
    }
    reads
}

pub fn dump_seq(id: &str, seq: &[u8], desc: Option<&str>) {
    let mut writer = fasta::Writer::new(io::stdout());
    writer.write(id, desc, seq).unwrap();
}

pub fn read2() {
    let s = b"ATCGATTCGATCGATTCGAT";
    println!("all {:?}", s);
    for w in s.windows(5) {
        println!("window {:?}", w);
    }
}
