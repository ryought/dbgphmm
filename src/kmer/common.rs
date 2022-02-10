//!
//! Kmer definitions
//!

pub trait NullableKmer {
    ///
    /// null <==> NNNNN
    ///
    fn is_null(&self) -> bool;
}

pub trait KmerLike:
    std::marker::Sized + PartialEq + NullableKmer + Eq + std::hash::Hash + Clone
{
    /// type of k+1-mer
    type Kp1mer: PartialEq + NullableKmer;
    /// type of k-1-mer
    type Km1mer: PartialEq + NullableKmer;
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
    fn prefix(&self) -> Self::Km1mer;
    ///
    /// suffix of the kmer
    /// ABBBB -> BBBB
    ///
    fn suffix(&self) -> Self::Km1mer;
    ///
    /// ABBBB and BBBBC is adjacent
    ///
    fn adjacent(&self, other: &Self) -> bool {
        self.suffix() == other.prefix()
    }
    ///
    /// XYYYY -> [YYYYA, YYYYC, YYYYG, YYYYT]
    ///
    fn childs(&self) -> Vec<Self>;
    ///
    /// YYYYZ -> [AYYYY, CYYYY, GYYYY, TYYYY]
    ///
    fn parents(&self) -> Vec<Self>;
    ///
    /// union of childs and parents
    /// XYYYZ -> [
    ///            YYYZA, YYYZC, YYYZG, YYYZT, (childs)
    ///            AXYYY, CXYYY, GXYYY, TXYYY  (parents)
    ///          ]
    ///
    fn neighbors(&self) -> Vec<Self> {
        let mut childs = self.childs();
        let parents = self.parents();
        childs.extend(parents);
        childs
    }
    ///
    /// XX (k mer) -> [AXX, CXX, GXX, TXX] (k+1 mer)
    ///
    fn preds(&self) -> Vec<Self::Kp1mer> {
        let bases = [b'A', b'C', b'G', b'T'];
        bases
            .iter()
            .map(|&first_base| self.extend_first(first_base))
            .collect()
    }
    ///
    /// XX (k mer) -> [XXA, XXC, XXG, XXT] (k+1 mer)
    ///
    fn succs(&self) -> Vec<Self::Kp1mer> {
        let bases = [b'A', b'C', b'G', b'T'];
        bases
            .iter()
            .map(|&last_base| self.extend_last(last_base))
            .collect()
    }
    ///
    /// head <==> NNNNX
    ///
    fn is_head(&self) -> bool {
        self.prefix().is_null()
    }
    ///
    /// tail <==> XNNNN
    ///
    fn is_tail(&self) -> bool {
        self.suffix().is_null()
    }
    ///
    /// last base is not N
    ///
    fn is_emitable(&self) -> bool {
        self.last() != b'N'
    }
    ///
    /// first base is N
    /// TODO check that the rest is not N
    ///
    fn is_starting(&self) -> bool {
        self.first() == b'N'
    }
    ///
    /// (YYY, X) -> XYYY
    ///
    fn extend_first(&self, first_base: u8) -> Self::Kp1mer;
    ///
    /// (YYY, Z) -> YYYZ
    ///
    fn extend_last(&self, last_base: u8) -> Self::Kp1mer;
    ///
    /// (XYYY, YYYZ) (two k mers) -> XYYYZ (k+1 mer)
    ///
    fn join(&self, other: &Self) -> Self::Kp1mer {
        if !self.adjacent(other) {
            panic!();
        }
        self.extend_last(other.last())
    }
    // construction
    // fn from(bases: &[u8]) -> Self;
    // fn to_vec(&self) -> Vec<u8>;
}

///
/// Most fundamental k-mer trait
/// TODO
///
pub trait KmerBase {
    fn k(&self) -> usize;
}

//
// Tests
//

#[cfg(test)]
mod tests {
    use super::*;
}
