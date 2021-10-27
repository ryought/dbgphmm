use crate::prob::Prob;
use itertools::izip;

#[derive(Debug, Clone)]
pub struct PHMMLayer {
    pub pM: Vec<Prob>,
    pub pI: Vec<Prob>,
    pub pD: Vec<Prob>,
    // Begin
    pub pMB: Prob,
    pub pIB: Prob,
    // End
    pub pE: Prob,
}

impl PHMMLayer {
    pub fn new(n_kmers: usize) -> PHMMLayer {
        PHMMLayer {
            pM: vec![Prob::from_prob(0.0); n_kmers],
            pI: vec![Prob::from_prob(0.0); n_kmers],
            pD: vec![Prob::from_prob(0.0); n_kmers],
            pMB: Prob::from_prob(0.0),
            pIB: Prob::from_prob(0.0),
            pE: Prob::from_prob(0.0),
        }
    }
    /// p[i][j] = (e[i] is emitted from j-th kmer)
    /// ignore 3 states of each kmer.
    pub fn to_kmer_prob(&self) -> Vec<Prob> {
        izip!(&self.pM, &self.pI, &self.pD)
            .map(|(&m, &i, &d)| m + i + d)
            .collect()
    }
    pub fn to_freqs(&self) -> Vec<f64> {
        self.to_kmer_prob()
            .into_iter()
            .map(|p| p.to_value())
            .collect()
    }
}

impl std::fmt::Display for PHMMLayer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Node:Begin\tpM={} pI={}", self.pMB, self.pIB);
        for i in 0..self.pD.len() {
            writeln!(
                f,
                "Node:{}\tpM={}\tpI={}\tpD={}",
                i, self.pM[i], self.pI[i], self.pD[i]
            );
        }
        writeln!(f, "Node:End\tpE={}", self.pE);
        Ok(())
    }
}

// operators
impl<'a, 'b> std::ops::Add<&'b PHMMLayer> for &'a PHMMLayer {
    type Output = PHMMLayer;
    fn add(self, other: &'b PHMMLayer) -> PHMMLayer {
        // assert that length is same
        assert_eq!(self.pM.len(), other.pM.len());
        assert_eq!(self.pI.len(), other.pI.len());
        assert_eq!(self.pD.len(), other.pD.len());
        let pM = self
            .pM
            .iter()
            .zip(other.pM.iter())
            .map(|(&s, &o)| s + o)
            .collect();
        let pI = self
            .pI
            .iter()
            .zip(other.pI.iter())
            .map(|(&s, &o)| s + o)
            .collect();
        let pD = self
            .pD
            .iter()
            .zip(other.pD.iter())
            .map(|(&s, &o)| s + o)
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

impl std::iter::Sum for PHMMLayer {
    fn sum<I: Iterator<Item = PHMMLayer>>(iter: I) -> PHMMLayer {
        iter.reduce(|a, b| &a + &b).unwrap()
    }
}
