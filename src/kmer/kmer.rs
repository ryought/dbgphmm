///
/// kmer base struct definitions
///

#[derive(Debug, PartialEq, PartialOrd, Eq, Hash, Clone)]
pub struct Kmer(Vec<u8>);

impl Kmer {
    pub fn from(s: &[u8]) -> Kmer {
        let v = s.to_vec();
        // assert items in v is a,c,g,t,n
        Kmer(v)
    }
    pub fn from_vec(v: Vec<u8>) -> Kmer {
        Kmer(v)
    }
    pub fn len(&self) -> usize {
        self.0.len()
    }
    pub fn adjacent(&self, other: &Kmer) -> bool {
        let (_, a_suffix) = self.0.split_first().expect("k should be >1");
        let (_, b_prefix) = other.0.split_last().expect("k should be >1");
        a_suffix == b_prefix
    }
    pub fn first(&self) -> u8 {
        let (first, _) = self.0.split_first().expect("k should be >=1");
        *first
    }
    pub fn last(&self) -> u8 {
        let (last, _) = self.0.split_last().expect("k should be >=1");
        *last
    }
    pub fn prefix(&self) -> Kmer {
        let (_, prefix) = self.0.split_last().expect("k should be >1");
        Kmer(prefix.to_vec())
    }
    pub fn suffix(&self) -> Kmer {
        let (_, suffix) = self.0.split_first().expect("k should be >1");
        Kmer(suffix.to_vec())
    }
    pub fn childs(&self) -> Vec<Kmer> {
        let (_, suffix) = self.0.split_first().unwrap();
        let childs = [b'A', b'C', b'G', b'T', b'N']
            .iter()
            .map(|last_base| {
                let mut v = suffix.to_vec();
                v.push(*last_base);
                Kmer::from_vec(v)
            })
            .collect();
        childs
    }
    pub fn parents(&self) -> Vec<Kmer> {
        let (_, prefix) = self.0.split_last().unwrap();
        let parents = [b'A', b'C', b'G', b'T', b'N']
            .iter()
            .map(|first_base| {
                let mut v = prefix.to_vec();
                v.insert(0, *first_base);
                Kmer::from_vec(v)
            })
            .collect();
        parents
    }
}

impl std::fmt::Display for Kmer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // iter returns reference
        for &b in self.0.iter() {
            write!(f, "{}", b as char)?;
        }
        Ok(())
    }
}

use std::collections::HashMap;
pub fn count(seq: &[u8], k: usize) -> HashMap<Kmer, usize> {
    let mut h: HashMap<Kmer, usize> = HashMap::new();
    for i in 0..=&seq.len() - k {
        let kmer = Kmer::from(&seq[i..i + k]);
        // let kmer = &seq[i..i + k];
        // let count = h.entry(kmer).or_insert(0);
        //*count += 1;
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

pub fn test() {
    let a = Kmer::from(b"ATCGATTAG");
    let b = Kmer::from(b"TCGATTAGT");
    let x = a.adjacent(&b);
    let y = a.last();
    println!("{} {} {} {} {}", a, b, a == b, x, y);
}

/// ATTC -> [NATT, NNAT, NNNA]
pub fn starting_kmers(kmer: &Kmer) -> Vec<Kmer> {
    let k = kmer.len();
    let blanks = std::iter::repeat(b'N').take(k).collect::<Vec<u8>>();
    let heads: Vec<Kmer> = (1..k)
        .rev()
        .map(|i| Kmer::from_vec([&blanks[i..], &kmer.0[..i]].concat()))
        .collect();
    heads
}

/// ATTC -> [TTCN, TCNN, CNNN]
pub fn ending_kmers(kmer: &Kmer) -> Vec<Kmer> {
    let k = kmer.len();
    let blanks = std::iter::repeat(b'N').take(k).collect::<Vec<u8>>();
    let tails: Vec<Kmer> = (1..k)
        .map(|i| Kmer::from_vec([&kmer.0[i..], &blanks[..i]].concat()))
        .collect();
    tails
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmer_equality() {
        let a = Kmer::from(b"ATCGATTAG");
        let b = Kmer::from(b"ATCGATTAG");
        assert_eq!(a, b);
    }
    #[test]
    fn kmer_adjacency() {
        let a = Kmer::from(b"ATCGATTAG");
        let b = Kmer::from(b"TCGATTAGA");
        assert!(a.adjacent(&b));
    }
    #[test]
    fn kmer_adjacency_fail() {
        let a = Kmer::from(b"ATCGATTAG");
        let b = Kmer::from(b"TCGATTAAA");
        assert!(!a.adjacent(&b));
    }
    #[test]
    fn kmer_last() {
        let a = Kmer::from(b"ATCGATTAG");
        assert_eq!(a.last(), b'G');
    }
    #[test]
    #[should_panic]
    fn k_zero_mer() {
        let a = Kmer::from(b"");
        let x = a.last();
    }
}
