//!
//! Histogram counter
//!
use crate::prob::Prob;
use fnv::FnvHashMap as HashMap;
use itertools::Itertools;

/// calculate (avarage, std dev, min, max) of the list of f64
pub fn stat(xs: &[f64]) -> (f64, f64, f64, f64) {
    let s: f64 = xs.iter().sum();
    let n: f64 = xs.len() as f64;
    if xs.len() == 0 {
        (0.0, 0.0, 0.0, 0.0)
    } else {
        let ave = s / n;
        let d: f64 = xs.iter().map(|x| (x - ave).powi(2)).sum();
        let sd = (d / n).sqrt();

        let min = xs
            .iter()
            .copied()
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        let max = xs
            .iter()
            .copied()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();

        (ave, sd, min, max)
    }
}

///
/// Represents (normalized) discrete probability distribution `P(X), X \in {0,1,2,...}`
///
#[derive(Clone, Debug)]
pub struct DiscreteDistribution {
    z: Prob,
    hashmap: HashMap<usize, Vec<Prob>>,
}

impl DiscreteDistribution {
    ///
    /// construct normalized distribution from unnormalized occurrences list
    /// convert a list of occurrences (x, p(x)) into {x: p(x)}
    ///
    pub fn from_occurs(probs: &[(usize, Prob)]) -> DiscreteDistribution {
        let mut h = HashMap::default();
        let mut z = Prob::zero();
        for &(x, p_x) in probs {
            h.entry(x).or_insert_with(|| Vec::new()).push(p_x);
            z += p_x;
        }
        DiscreteDistribution { hashmap: h, z }
    }
    pub fn to_normalized(&self) -> Vec<(usize, Prob)> {
        self.hashmap.keys().map(|&x| (x, self.p_x(x))).collect()
    }
    ///
    /// Get P(X=x)
    ///
    pub fn p_x(&self, x: usize) -> Prob {
        match self.hashmap.get(&x) {
            Some(probs) => {
                let p: Prob = probs.iter().sum();
                p / self.z
            }
            None => Prob::zero(),
        }
    }
    ///
    /// E[X] = sum_x(x * P(X=x))
    ///
    pub fn mean(&self) -> f64 {
        self.hashmap
            .keys()
            .map(|&x| x as f64 * self.p_x(x).to_value())
            .sum()
    }
    pub fn to_short_string(&self) -> String {
        self.hashmap
            .keys()
            .sorted()
            .map(|&x| format!("p({})={:.3}", x, self.p_x(x).to_value()))
            .join(",")
    }
}
impl std::fmt::Display for DiscreteDistribution {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.hashmap
                .keys()
                .sorted()
                .map(|&x| format!("p(x={})={:.5}", x, self.p_x(x)))
                .join(",")
        )
    }
}

///
/// Histogram counter struct
///
/// * new
/// * add
/// * get
/// * range
/// * iter
///
#[derive(Clone, Debug)]
pub struct Hist(HashMap<usize, usize>);

impl Hist {
    ///
    /// Create an empty histogram counter
    ///
    pub fn new() -> Self {
        Hist(HashMap::default())
    }
    ///
    /// create histogram from occurance table
    ///
    pub fn from(values: &[usize]) -> Self {
        let mut h = Hist::new();
        for &value in values {
            h.add(value);
        }
        h
    }
    ///
    /// Increment a count of the value
    ///
    pub fn add(&mut self, value: usize) {
        *self.0.entry(value).or_insert(0) += 1;
    }
    ///
    /// Get count of the value
    ///
    pub fn get(&self, value: usize) -> usize {
        match self.0.get(&value) {
            Some(&count) => count,
            None => 0,
        }
    }
    ///
    /// Get the min/max of the values
    ///
    pub fn range(&self) -> Option<(usize, usize)> {
        if !self.0.is_empty() {
            let min = *self.0.keys().min().unwrap();
            let max = *self.0.keys().max().unwrap();
            Some((min, max))
        } else {
            None
        }
    }
    ///
    /// Iterate over values and its counts
    ///
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.0.iter().map(|(k, v)| (*k, *v))
    }
    ///
    /// Get the number of elements stored in the histogram
    ///
    pub fn len(&self) -> usize {
        self.0.values().sum()
    }
    ///
    /// the histgram is empty or not.
    ///
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl std::fmt::Display for Hist {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.iter()
            .sorted()
            .enumerate()
            .try_for_each(|(i, (k, v))| write!(f, "{}{}:{}", if i != 0 { "," } else { "" }, k, v))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prob::p;
    #[test]
    fn hist_test() {
        let mut h = Hist::new();
        h.add(10);
        h.add(10);
        h.add(2);
        h.add(10);
        h.add(2);
        h.add(3);
        println!("{:?}", h);
        assert_eq!(h.range(), Some((2, 10)));
        assert_eq!(h.get(10), 3);
        assert_eq!(h.get(2), 2);
        assert_eq!(h.get(3), 1);
        assert_eq!(h.get(0), 0);
        let c: Vec<_> = h.iter().collect();
        assert_eq!(c, vec![(2, 2), (3, 1), (10, 3)]);
        println!("{}", h);
        assert_eq!(h.to_string(), "2:2,3:1,10:3");
    }
    #[test]
    fn hist_test_from() {
        let h = Hist::from(&[9, 8, 2, 2, 8, 5]);
        println!("{}", h);
        assert_eq!(h.to_string(), "2:2,5:1,8:2,9:1");
    }
    #[test]
    fn stat_ave_sd() {
        let (ave, sd, min, max) = stat(&vec![1., 1., 1., 1., 1., 1., 1.]);
        assert_eq!(ave, 1.0);
        assert_eq!(sd, 0.0);
        assert_eq!(min, 1.0);
        assert_eq!(max, 1.0);
        println!("{} {} {} {}", ave, sd, min, max);

        let (ave, sd, min, max) = stat(&vec![0., 10.]);
        assert_eq!(ave, 5.0);
        assert_eq!(sd, 5.0);
        assert_eq!(min, 0.0);
        assert_eq!(max, 10.0);
        println!("{} {} {} {}", ave, sd, min, max);
    }
    #[test]
    fn hist_get_normalized() {
        let d = DiscreteDistribution::from_occurs(&[
            (0, p(0.5)),
            (0, p(0.2)),
            (1, p(0.2)),
            (0, p(0.1)),
        ]);
        assert_abs_diff_eq!(d.p_x(0), p(0.8));
        assert_abs_diff_eq!(d.p_x(1), p(0.2));
        println!("E[X]={}", d.mean());
        assert_abs_diff_eq!(d.mean(), 0.2);

        let d = DiscreteDistribution::from_occurs(&[(0, p(0.5)), (0, p(0.2))]);
        assert_abs_diff_eq!(d.p_x(0), p(1.0));
        println!("E[X]={}", d.mean());
        assert_abs_diff_eq!(d.mean(), 0.0);

        let d = DiscreteDistribution::from_occurs(&[(0, p(0.3)), (1, p(0.1)), (1, p(0.1))]);
        assert_abs_diff_eq!(d.p_x(0), p(0.6));
        assert_abs_diff_eq!(d.p_x(1), p(0.4), epsilon = 0.000001);
        println!("E[X]={}", d.mean());
        assert_abs_diff_eq!(d.mean(), 0.4);
    }
}
