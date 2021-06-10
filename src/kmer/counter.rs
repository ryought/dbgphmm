///
/// kmer-counters
///
extern crate rand;
use bio::io::fasta;
use rand::seq::SliceRandom;
use rand::{thread_rng, Rng};

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

pub fn test_counter() -> usize {
    let s = get_random_vec(1000000);
    // println!("random {}", print(&s));
    // let h = count(b"ATGCTAGCTTATG", 3);
    let h = count(&s, 8);
    println!("{}", h.len());
    /*
    for (i, (kmer, &occ)) in h.iter().enumerate() {
        println!("{}: {} {}", i, print(kmer), occ);
    }
    */
    h.len()
}

pub fn get_random_vec(len: usize) -> Vec<u8> {
    // let mut rng = rand::thread_rng();
    // let v: Vec<u8> = bases.choose_multiple(&mut rng, len).cloned().collect();
    let v: Vec<u8> = Vec::new();
    // let choices = [1, 2, 4, 8, 16, 32];
    let choices = "ACGT".as_bytes();
    let mut rng = rand::thread_rng();
    // println!("{:?}", choices.choose(&mut rng));
    let v: Vec<u8> = std::iter::repeat(())
        .map(|()| *choices.choose(&mut rng).unwrap())
        .take(len)
        .collect();
    v
}

use std::collections::HashMap;
pub fn count(seq: &[u8], k: usize) -> HashMap<Vec<u8>, usize> {
    let mut h = HashMap::new();
    for i in 0..=&seq.len() - k {
        let kmer = &seq[i..i + k];
        match h.get(kmer) {
            Some(&occ) => {
                h.insert(kmer.to_vec(), occ + 1);
            }
            _ => {
                h.insert(kmer.to_vec(), 1);
            }
        };
    }
    h
}

use super::kmer::Kmer;
pub fn count_my_kmer(seq: &[u8], k: usize) -> HashMap<Kmer, usize> {
    let mut h = HashMap::new();
    for i in 0..=&seq.len() - k {
        // XXX bottleneck here
        let kmer = Kmer::from(&seq[i..i + k]);
        match h.get(&kmer) {
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
