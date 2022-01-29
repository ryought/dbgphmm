//!
//! Kmer definitions
//!
use super::kmer::Kmer;

pub trait KmerLike: std::marker::Sized {
    ///
    /// k of the k-mer
    ///
    fn len(&self) -> usize;
    ///
    /// k of the k-mer
    /// (an alias of KmerLike.len)
    ///
    fn k(&self) -> usize {
        self.len()
    }
    ///
    /// ABBBB -> A
    ///
    fn first(&self) -> u8;
    ///
    /// AAAAB -> B
    ///
    fn last(&self) -> u8;
    ///
    /// prefix of the kmer
    /// AAAAB -> AAAA
    ///
    fn prefix(&self) -> Kmer;
    ///
    /// suffix of the kmer
    /// ABBBB -> BBBB
    ///
    fn suffix(&self) -> Kmer;
    ///
    /// ABBBB and BBBBC is adjacent
    ///
    fn adjacent(&self, other: &Kmer) -> bool;
    // TODO default implementation
    // fn adjacent(&self, other: &Kmer) -> bool {
    //     self.suffix() == other.prefix()
    // }
    ///
    /// XYYYY -> [YYYYA, YYYYC, YYYYG, YYYYT]
    ///
    fn childs(&self) -> Vec<Kmer>;
    ///
    /// YYYYZ -> [AYYYY, CYYYY, GYYYY, TYYYY]
    ///
    fn parents(&self) -> Vec<Kmer>;
    ///
    /// union of childs and parents
    /// XYYYZ -> [
    ///            YYYZA, YYYZC, YYYZG, YYYZT, (childs)
    ///            AXYYY, CXYYY, GXYYY, TXYYY  (parents)
    ///          ]
    ///
    fn neighbors(&self) -> Vec<Kmer>;
    ///
    /// XX (k mer) -> [AXX, CXX, GXX, TXX] (k+1 mer)
    ///
    fn preds(&self) -> Vec<Kmer>;
    ///
    /// XX (k mer) -> [XXA, XXC, XXG, XXT] (k+1 mer)
    ///
    fn succs(&self) -> Vec<Kmer>;
    ///
    /// (XYYY, YYYZ) (two k mers) -> XYYYZ (k+1 mer)
    ///
    fn join(&self, other: &Kmer) -> Kmer;
    ///
    /// head <==> NNNNX
    ///
    fn is_head(&self) -> bool;
    ///
    /// tail <==> XNNNN
    ///
    fn is_tail(&self) -> bool;
    ///
    /// last base is not N
    ///
    fn is_emitable(&self) -> bool;
    ///
    /// first base is N
    ///
    fn is_starting(&self) -> bool;
    // internal functions
    fn extend_first(&self, first_base: u8) -> Kmer;
    fn extend_last(&self, last_base: u8) -> Kmer;
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
