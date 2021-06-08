///
/// kmer-counters
///
extern crate bio;
use bio::io::fasta;

pub struct Config {
    pub filename: String,
}
fn parse_config(args: &[String]) -> Config {
    Config {
        filename: args[1].clone(),
    }
}
impl Config {
    pub fn new(args: &[String]) -> Config {
        if args.len() < 2 {
            panic!("not enough arguments");
        }
        Config {
            filename: args[1].clone(),
        }
    }
}

fn print(kmer: &[u8]) -> &str {
    std::str::from_utf8(kmer).unwrap()
}

pub fn run_counter(config: Config) {
    println!("filename: {}", config.filename);
    let mut reader = fasta::Reader::from_file(config.filename).unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let h = count(record.seq(), 20);
        for (kmer, &occ) in h.iter() {
            println!("{} {}", print(kmer), occ);
        }
    }
}

pub fn test_counter() {
    let h = count(b"ATGCTAGCTTATG", 3);
    for (kmer, &occ) in h.iter() {
        println!("{} {}", print(kmer), occ);
    }
}

use std::collections::HashMap;
fn count(seq: &[u8], k: usize) -> HashMap<&[u8], usize> {
    let mut h = HashMap::new();
    for i in 0..=&seq.len() - k {
        let kmer = &seq[i..i + k];
        match h.get(kmer) {
            Some(&occ) => {
                h.insert(kmer, occ + 1);
            }
            _ => {
                h.insert(kmer, 1);
            }
        };
    }
    h
}
