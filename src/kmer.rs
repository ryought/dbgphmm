///
/// kmer base struct definitions
///

#[derive(Debug, PartialEq, PartialOrd, Eq, Hash)]
pub struct Kmer(Vec<u8>);

impl Kmer {
    pub fn from(s: &[u8]) -> Kmer {
        let v = s.to_vec();
        // assert items in v is a,c,g,t,n
        Kmer(v)
    }
    pub fn adjacent(&self, other: &Kmer) -> bool {
        let (_, a_suffix) = self.0.split_first().expect("k should be >1");
        let (_, b_prefix) = other.0.split_last().expect("k should be >1");
        a_suffix == b_prefix
    }
    pub fn last(&self) -> u8 {
        let (last, _) = self.0.split_last().expect("k should be >=1");
        *last
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
