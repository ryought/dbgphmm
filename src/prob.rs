///
/// probability calculation
/// implements logaddexp
///
use approx::AbsDiffEq;
use serde_with::{DeserializeFromStr, SerializeDisplay};
use std::str::FromStr;

///
/// Wrapper of f64 that represents probability `0 <= p <= 1`
///
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, SerializeDisplay, DeserializeFromStr)]
pub struct Prob(f64);

///
/// short-hand of `Prob::from_prob`
///
pub fn p(p: f64) -> Prob {
    Prob::from_prob(p)
}

///
/// short-hand of `Prob::from_log_prob`
///
pub fn lp(lp: f64) -> Prob {
    Prob::from_log_prob(lp)
}

use once_cell::sync::Lazy;
const MAX_PRECALCULATED_X: usize = 10;
static ys: Lazy<[f64; MAX_PRECALCULATED_X]> = Lazy::new(|| {
    let mut v = [0f64; MAX_PRECALCULATED_X];
    for x in 0..MAX_PRECALCULATED_X {
        v[x] = (x as f64).ln();
    }
    v
});

///
/// Faster cached log function `(x as f64).ln()`
///
#[inline]
pub fn ln_int(x: usize) -> f64 {
    // if x is small, return the precalculated
    if x < MAX_PRECALCULATED_X {
        ys[x]
    } else {
        (x as f64).ln()
    }
}

pub fn ln_int0(x: usize) -> f64 {
    (x as f64).ln()
}

impl Prob {
    ///
    ///
    pub fn from_prob(value: f64) -> Prob {
        Prob(value.ln())
    }
    ///
    ///
    pub fn from_log_prob(log_value: f64) -> Prob {
        Prob(log_value)
    }
    ///
    /// Get the probability (in `[0, 1]`)
    pub fn to_value(self) -> f64 {
        self.0.exp()
    }
    ///
    /// Get the log probability
    pub fn to_log_value(self) -> f64 {
        self.0
    }
    ///
    /// Is `p == 0` or not? (log p = -inf)
    ///
    pub fn is_zero(self) -> bool {
        self.0.is_infinite() && self.0.is_sign_negative()
    }
    ///
    /// Is `p == 1`? (log p = 0)
    ///
    pub fn is_one(self) -> bool {
        self.0 == 0.0
    }
    ///
    /// prob=0.0
    ///
    pub fn zero() -> Prob {
        Prob(f64::NEG_INFINITY)
    }
    ///
    /// prob=1.0
    ///
    pub fn one() -> Prob {
        Prob(0.0)
    }
    ///
    /// abs diff of two probs `= |p_a - p_b|`
    ///
    pub fn diff(&self, other: Prob) -> f64 {
        (self.to_value() - other.to_value()).abs()
    }
    ///
    /// abs diff of two log probs `= |log p_a - log p_b|`
    ///
    pub fn log_diff(&self, other: Prob) -> f64 {
        if self.is_zero() {
            if other.is_zero() {
                0.0
            } else {
                f64::INFINITY
            }
        } else {
            if other.is_zero() {
                f64::INFINITY
            } else {
                (self.to_log_value() - other.to_log_value()).abs()
            }
        }
    }
}

/// p=0 (Prob(-inf)) as a default value
impl Default for Prob {
    fn default() -> Self {
        Prob(f64::NEG_INFINITY)
    }
}

///
/// Prob has multiplicative identity element
/// `num_traits::One`
///
impl num_traits::One for Prob {
    fn one() -> Self {
        Prob::one()
    }
}

///
/// Prob has additive identity element
/// `num_traits::Zero`
///
impl num_traits::Zero for Prob {
    fn zero() -> Self {
        Prob::zero()
    }
    fn is_zero(&self) -> bool {
        Prob::is_zero(*self)
    }
}

// display
impl std::fmt::Display for Prob {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}({:.4})", self.0, self.to_value())
    }
}
impl FromStr for Prob {
    type Err = std::num::ParseFloatError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (front, _) = s.split_once('(').unwrap();
        front.parse::<f64>().map(|p| Prob(p))
    }
}

/// Addition of two probabilities `px + py` in log space
///
/// If `px > py`:
///
/// ```text
/// log(exp(x) + exp(y))
///  = log(exp(x) (1 + exp(y-x)))
///  = log(exp(x)) + log(1 + exp(y-x))
///  = x + log(1 + exp(y-x))
/// ```
impl std::ops::Add for Prob {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        let x = self.0;
        let y = other.0;
        let (x, y) = if x >= y { (x, y) } else { (y, x) };
        if y == f64::NEG_INFINITY {
            // x + 0 = x
            Prob(x)
        } else if x == y {
            // x + x = 2x
            Prob(x + 2f64.ln())
        } else {
            Prob(x + (y - x).exp().ln_1p())
        }
    }
}

/// Multiplication of two probabilities `px * py` in log space
///
/// ```text
/// log(px * py) = log(px) + log(py)
/// ```
impl std::ops::Mul for Prob {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Prob(self.0 + other.0)
    }
}

/// Division of two probabilities `px / py` in log space
///
/// ```text
/// log(px / py) = log(px) - log(py)
/// ```
impl std::ops::Div for Prob {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Prob(self.0 - other.0)
    }
}

// assign
impl std::ops::AddAssign for Prob {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}
impl std::ops::MulAssign for Prob {
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}
// sum/prod
impl std::iter::Sum for Prob {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Prob::from_prob(0.0), |a, b| a + b)
    }
}
impl<'a> std::iter::Sum<&'a Self> for Prob {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Prob::from_prob(0.0), |a, b| a + *b)
    }
}
impl std::iter::Product for Prob {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Prob::from_prob(1.0), |a, b| a * b)
    }
}
impl<'a> std::iter::Product<&'a Self> for Prob {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Prob::from_prob(1.0), |a, b| a * *b)
    }
}

//
// Prob mul/div usize
//

/// Multiplication of Prob and usize `p * c`
///
impl std::ops::Mul<usize> for Prob {
    type Output = Self;
    fn mul(self, rhs: usize) -> Self {
        Prob(self.0 + (rhs as f64).ln())
    }
}

/// Division of Prob and usize `p / c`
///
impl std::ops::Div<usize> for Prob {
    type Output = Self;
    fn div(self, rhs: usize) -> Self {
        if rhs == 0 {
            panic!("zero division error")
        } else {
            Prob(self.0 - (rhs as f64).ln())
        }
    }
}

/// for approx `assert_abs_diff_eq`
impl AbsDiffEq for Prob {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        f64::abs_diff_eq(&self.0, &other.0, epsilon)
    }
}

impl Eq for Prob {}
impl Ord for Prob {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_id() {
        let x = Prob::from_prob(0.3);
        let e = Prob::from_prob(0.0);
        assert_relative_eq!((x + e).0, x.0);
        assert_relative_eq!((x * e).0, e.0);
    }
    #[test]
    fn test_id_id() {
        let x = Prob::from_prob(1.0);
        let e = Prob::from_prob(0.0);
        assert_relative_eq!((x + e).0, x.0);
        assert_relative_eq!((x * e).0, e.0);
    }
    #[test]
    fn test_sum() {
        let xs = vec![
            Prob::from_prob(0.1),
            Prob::from_prob(0.1),
            Prob::from_prob(0.1),
        ];
        let x: Prob = xs.iter().sum();
        let y = Prob::from_prob(0.3);
        assert_relative_eq!(x.to_value(), y.to_value());
    }
    #[test]
    fn test_prod() {
        let xs = vec![
            Prob::from_prob(0.1),
            Prob::from_prob(0.1),
            Prob::from_prob(0.1),
        ];
        let x: Prob = xs.iter().product();
        let y = Prob::from_prob(0.001);
        assert_relative_eq!(x.to_value(), y.to_value());
    }
    #[test]
    fn prob_add_mul() {
        assert_eq!(p(0.0) + p(1.0), p(1.0));
        assert_eq!(p(0.0) * p(1.0), p(0.0));
        assert_abs_diff_eq!((p(0.3) + p(0.3)).0, p(0.6).0);
        assert_abs_diff_eq!((p(0.3) * p(0.3)).0, p(0.09).0);
        assert_abs_diff_eq!((p(0.5) + p(0.00001)).0, p(0.50001).0);
        assert_abs_diff_eq!((p(0.5) * p(0.00001)).0, p(0.000005).0);
    }
    #[test]
    fn prob_sum_prod() {
        // sum/prod of zero element vec
        let xs: Vec<Prob> = vec![];
        let sum: Prob = xs.iter().sum();
        let product: Prob = xs.iter().product();
        assert_eq!(sum, p(0.0));
        assert_eq!(product, p(1.0));

        // sum/prod of vec of p=0
        let xs: Vec<Prob> = vec![p(0.0), p(0.0)];
        let sum: Prob = xs.iter().sum();
        let product: Prob = xs.iter().product();
        assert_eq!(sum, p(0.0));
        assert_eq!(product, p(0.0));
    }
    #[test]
    fn test_reflect() {
        let x = Prob::from_prob(1.0);
        let e = Prob::from_prob(0.0);
        assert_relative_eq!((x + e).0, x.0);
        assert_relative_eq!((x * e).0, e.0);
        assert_relative_eq!((e * e).0, e.0);
    }
    #[test]
    fn test_zero() {
        let zero = Prob::from_prob(0.0);
        println!("{:?}", zero);
        assert!(zero.is_zero());
        let nonzero = Prob::from_prob(0.00001);
        assert!(!nonzero.is_zero());
    }
    #[test]
    fn test_prob_assign() {
        let mut x = p(0.4);
        let y = p(0.2);
        x += y;
        assert_abs_diff_eq!(x, p(0.6));
        let z = p(0.5);
        x *= z;
        assert_abs_diff_eq!(x, p(0.3));
        let o = p(1.0);
        x *= o;
        assert_abs_diff_eq!(x, p(0.3));
        let z = p(0.0);
        x += z;
        assert_abs_diff_eq!(x, p(0.3));
        x *= z;
        assert!(x.is_zero());
    }
    #[test]
    fn prob_sort() {
        // Sort by Ord and Eq
        let mut ps = vec![p(0.9), p(0.2), p(0.5), p(0.1), p(1.0), p(0.0)];
        ps.sort();
        println!("{:?}", ps);
        assert_eq!(ps[0], p(0.0));
        assert_eq!(ps[1], p(0.1));
        assert_eq!(ps[2], p(0.2));
        assert_eq!(ps[3], p(0.5));
        assert_eq!(ps[4], p(0.9));
        assert_eq!(ps[5], p(1.0));
    }
    #[test]
    fn prob_max_min() {
        let mut ps = vec![p(0.9), p(0.2), p(0.5), p(0.1), p(1.0), p(0.0)];
        let max = ps.iter().max().unwrap();
        assert_eq!(*max, p(1.0));
        let min = ps.iter().min().unwrap();
        assert_eq!(*min, p(0.0));

        assert!(p(0.1) > p(0.09999));
        assert!(p(0.1) < p(0.100001));
        assert!(p(0.0) < p(0.01));
        assert!(p(1.0) > p(0.01));
    }
    #[test]
    fn prob_assert_eq() {
        assert!(abs_diff_eq!(p(0.1), p(0.1)));
        assert!(!abs_diff_eq!(p(0.1), p(0.2)));
        assert!(!abs_diff_eq!(p(0.1), p(0.11)));
        assert!(abs_diff_eq!(p(0.1), p(0.11), epsilon = 0.1));
        assert!(abs_diff_eq!(p(1.0), p(1.0)));
        // TODO
        assert!(!abs_diff_eq!(p(0.0), p(0.0)));
    }
    #[test]
    fn prob_zero_one() {
        assert_eq!(Prob::one(), Prob::from_prob(1.0));
        assert_eq!(Prob::zero(), Prob::from_prob(0.0));
        assert!(Prob::zero().is_zero());
    }
    #[test]
    fn prob_serialize() {
        // Display and FromStr
        let p1 = Prob::one();
        let p05 = Prob::from_prob(0.5);
        let p0 = Prob::zero();
        println!("{} {} {}", p1, p05, p0);
        assert_eq!(Prob::from_str(&p1.to_string()).unwrap(), p1);
        assert_eq!(Prob::from_str(&p05.to_string()).unwrap(), p05);
        assert_eq!(Prob::from_str(&p0.to_string()).unwrap(), p0);

        let f = |p: Prob| {
            let json = &serde_json::to_string(&p).unwrap();
            println!("p={} json={}", p, json);
            serde_json::from_str(&json).unwrap()
        };
        assert_eq!(p1, f(p1));
        assert_eq!(p05, f(p05));
        assert_eq!(p0, f(p0));
    }
    #[test]
    fn prob_diff() {
        let p1 = Prob::one();
        let p05 = Prob::from_prob(0.5);
        let p06 = Prob::from_prob(0.6);
        let p0 = Prob::zero();
        assert_eq!(0.0, p1.log_diff(p1));
        assert_eq!(0.0, p05.log_diff(p05));
        assert_eq!(0.0, p06.log_diff(p06));
        assert_eq!(0.0, p0.log_diff(p0));

        assert_eq!(f64::INFINITY, p0.log_diff(p1));
        assert_eq!(f64::INFINITY, p1.log_diff(p0));
    }
    #[test]
    fn prob_muldiv_usize() {
        assert_eq!(p(0.0), p(1.0) * 0);
        assert_eq!(p(1.0), p(1.0) * 1);
        assert_eq!(p(2.0), p(1.0) * 2);
        assert_eq!(p(1.0), p(1.0) / 1);
        assert_eq!(p(0.5), p(1.0) / 2);

        assert_eq!(p(0.0), p(0.5) * 0);
        assert_eq!(p(0.5), p(0.5) * 1);
        assert_eq!(p(1.0), p(0.5) * 2);
        assert_eq!(p(0.5), p(0.5) / 1);
        assert_eq!(p(0.1), p(0.5) / 5);

        assert_eq!(p(0.0), p(0.0) * 0);
        assert_eq!(p(0.0), p(0.0) * 1);
        assert_eq!(p(0.0), p(0.0) * 2);
        assert_eq!(p(0.0), p(0.0) / 1);
        assert_eq!(p(0.0), p(0.0) / 2);
    }
    #[test]
    fn const_log_int() {
        for x in 0..100 {
            println!("{}", ln_int(x));
            assert_eq!(ln_int(x), (x as f64).ln());
        }
    }
    use test::Bencher;
    #[bench]
    fn usize_log_cached(b: &mut Bencher) {
        // precalculation
        test::black_box(ln_int(0));
        b.iter(|| {
            let mut r = 0.0;
            for x in 0..10 {
                // r += ln_int(test::black_box(x));
                r += ys[x];
            }
        })
    }
    #[bench]
    fn usize_log_normal(b: &mut Bencher) {
        b.iter(|| {
            let mut r = 0.0;
            for x in 0..10 {
                r += (test::black_box(x) as f64).ln();
            }
        })
    }

    use rand::prelude::*;
    use rand_xoshiro::Xoshiro256PlusPlus;
    #[bench]
    fn prob_add_zero(b: &mut Bencher) {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
        let n = 100;
        let m = 10;
        let mut xs = vec![Prob::zero(); n];
        for _ in 0..m {
            let i: usize = rng.gen_range(1..n);
            xs[i] = Prob::one();
        }
        println!("xs={:?}", xs);
        let mut r = Prob::one();
        b.iter(|| {
            for i in 1..n {
                r += xs[i];
            }
        });
        println!("r={}", r);
    }
}
