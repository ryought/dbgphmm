//!
//! Kmer definitions
//!
use super::kmer::Kmer;

pub trait KmerLike: std::marker::Sized {
    fn len(&self) -> usize;
    fn k(&self) -> usize {
        self.len()
    }
    fn adjacent(&self, other: &Kmer) -> bool;
    fn first(&self) -> u8;
    fn last(&self) -> u8;
    fn prefix(&self) -> Kmer;
    fn suffix(&self) -> Kmer;
    fn childs(&self) -> Vec<Kmer>;
    fn parents(&self) -> Vec<Kmer>;
    fn neighbors(&self) -> Vec<Kmer>;
    fn preds(&self) -> Vec<Kmer>;
    fn succs(&self) -> Vec<Kmer>;
    fn extend_first(&self, first_base: u8) -> Kmer;
    fn extend_last(&self, last_base: u8) -> Kmer;
    fn join(&self, other: &Kmer) -> Kmer;
    fn is_head(&self) -> bool;
    fn is_tail(&self) -> bool;
    fn is_emitable(&self) -> bool;
    fn is_starting(&self) -> bool;
    // construction
    // fn from(bases: &[u8]) -> Self;
    // fn to_vec(&self) -> Vec<u8>;
}

///
/// Most fundamental k-mer trait
///
pub trait KmerBase {
    type Kp1mer;
    type Km1mer;
    fn k(&self) -> usize;
    fn prefix(&self) -> Self::Km1mer;
    fn suffix(&self) -> Self::Km1mer;
}

//
// Tests
//

#[cfg(test)]
mod tests {
    use super::*;
}
