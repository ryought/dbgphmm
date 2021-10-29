use crate::prob::Prob;
use crate::veclike::{DenseVec, VecLike};
use itertools::izip;

#[derive(Debug, Clone)]
pub struct PHMMLayer<V: VecLike<Prob>> {
    pub pM: V,
    pub pI: V,
    pub pD: V,
    // Begin
    pub pMB: Prob,
    pub pIB: Prob,
    // End
    pub pE: Prob,
}

impl<V: VecLike<Prob>> PHMMLayer<V> {
    pub fn new(n_kmers: usize) -> Self {
        PHMMLayer {
            pM: V::new(n_kmers, Prob::from_prob(0.0)),
            pI: V::new(n_kmers, Prob::from_prob(0.0)),
            pD: V::new(n_kmers, Prob::from_prob(0.0)),
            pMB: Prob::from_prob(0.0),
            pIB: Prob::from_prob(0.0),
            pE: Prob::from_prob(0.0),
        }
    }
    /// p[i][j] = (e[i] is emitted from j-th kmer)
    /// ignore 3 states of each kmer.
    pub fn to_kmer_prob(&self) -> Vec<Prob> {
        izip!(self.pM.iter(), self.pI.iter(), self.pD.iter())
            .map(|(m, i, d)| m + i + d)
            .collect()
    }
    pub fn to_freqs(&self) -> Vec<f64> {
        self.to_kmer_prob()
            .into_iter()
            .map(|p| p.to_value())
            .collect()
    }
}

impl<V: VecLike<Prob>> std::fmt::Display for PHMMLayer<V> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Node:Begin\tpM={} pI={}", self.pMB, self.pIB);
        for i in 0..self.pD.len() {
            writeln!(
                f,
                "Node:{}\tpM={}\tpI={}\tpD={}",
                i,
                self.pM.get(i),
                self.pI.get(i),
                self.pD.get(i)
            );
        }
        writeln!(f, "Node:End\tpE={}", self.pE);
        Ok(())
    }
}

// operators
impl<'a, 'b> std::ops::Add<&'b PHMMLayer<DenseVec<Prob>>> for &'a PHMMLayer<DenseVec<Prob>> {
    type Output = PHMMLayer<DenseVec<Prob>>;
    fn add(self, other: &'b PHMMLayer<DenseVec<Prob>>) -> PHMMLayer<DenseVec<Prob>> {
        // assert that length is same
        assert_eq!(self.pM.len(), other.pM.len());
        assert_eq!(self.pI.len(), other.pI.len());
        assert_eq!(self.pD.len(), other.pD.len());
        let pM = self
            .pM
            .iter()
            .zip(other.pM.iter())
            .map(|(s, o)| s + o)
            .collect();
        let pI = self
            .pI
            .iter()
            .zip(other.pI.iter())
            .map(|(s, o)| s + o)
            .collect();
        let pD = self
            .pD
            .iter()
            .zip(other.pD.iter())
            .map(|(s, o)| s + o)
            .collect();
        PHMMLayer {
            pM,
            pI,
            pD,
            pMB: self.pMB + other.pMB,
            pIB: self.pIB + other.pIB,
            pE: self.pE + other.pE,
        }
    }
}

impl std::iter::Sum for PHMMLayer<DenseVec<Prob>> {
    fn sum<I: Iterator<Item = PHMMLayer<DenseVec<Prob>>>>(iter: I) -> PHMMLayer<DenseVec<Prob>> {
        iter.reduce(|a, b| &a + &b).unwrap()
    }
}
