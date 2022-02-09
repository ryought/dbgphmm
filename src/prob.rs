///
/// probability calculation
/// implements logaddexp
///
use approx::AbsDiffEq;

///
/// Wrapper of f64 that represents probability `0 <= p <= 1`
///
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
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

impl Prob {
    pub fn from_prob(value: f64) -> Prob {
        Prob(value.ln())
    }
    pub fn from_log_prob(log_value: f64) -> Prob {
        Prob(log_value)
    }
    pub fn to_value(self) -> f64 {
        self.0.exp()
    }
    pub fn to_log_value(self) -> f64 {
        self.0
    }
    ///
    /// Is `p == 0` or not?
    pub fn is_zero(self) -> bool {
        self.0.is_infinite() && self.0.is_sign_negative()
    }
}

// display
impl std::fmt::Display for Prob {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:.2}(=log({}))", self.0, self.to_value())
    }
}

// operators
impl std::ops::Add for Prob {
    type Output = Self;
    /// Given x>y
    /// log(exp(x) + exp(y))
    /// = log(exp(x) (1 + exp(y-x)))
    /// = log(exp(x)) + log(1 + exp(y-x))
    /// = x + log(1 + exp(y-x))
    fn add(self, other: Self) -> Self {
        let x = self.0;
        let y = other.0;
        if x == y {
            Prob(x + 2f64.ln())
        } else if x > y {
            Prob(x + ((y - x).exp() + 1.0).ln())
        } else {
            Prob(y + ((x - y).exp() + 1.0).ln())
        }
    }
}
impl std::ops::Mul for Prob {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Prob(self.0 + other.0)
    }
}
impl std::ops::Div for Prob {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Prob(self.0 - other.0)
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

pub fn test() {
    let z0 = Prob::from_prob(0.0);
    let z1 = Prob::from_prob(1.0);
    let x = Prob::from_prob(0.3);
    let y = Prob::from_prob(0.3);
    let _a2 = Prob::from_prob(0.6);
    let a = x + y;
    let b2 = Prob::from_prob(0.09);
    let b = x * y;
    println!(
        "{} {} {} {} {} {} {}",
        x.to_value(),
        y.to_value(),
        a.to_value(),
        b.to_value(),
        z0.0,
        z1.0,
        (z0 + z1).to_value(),
    );
    assert_abs_diff_eq!(b2.0, b.0);

    let x = Prob::from_prob(0.3);
    let e = Prob::from_prob(0.0);
    assert_relative_eq!((x + e).0, x.0);
    assert_relative_eq!((x * e).0, e.0);
    println!("{} == {}", (x + e).to_value(), x.to_value());
    println!("{} == {}", (x * e).to_value(), e.to_value());

    // sum
    let xs = vec![
        Prob::from_prob(0.1),
        Prob::from_prob(0.1),
        Prob::from_prob(0.1),
        Prob::from_prob(0.1),
    ];
    let s: Prob = xs.iter().sum();
    let p: Prob = xs.iter().product();
    let x = xs.iter().fold(Prob::from_prob(0.0), |sum, i| sum + *i);
    let y = Prob::from_prob(0.4);
    assert_relative_eq!(x.0, y.0);
    println!(
        "{} {} {} {}",
        x.to_value(),
        s.to_value(),
        y.to_value(),
        p.to_value()
    );

    let xs: Vec<Prob> = vec![];
    let s: Prob = xs.iter().sum();
    let p: Prob = xs.iter().product();
    println!("0={} 1={}", s, p);

    let xs: Vec<Prob> = vec![Prob::from_prob(0.0)];
    let s: Prob = xs.iter().sum();
    let p: Prob = xs.iter().product();
    println!("0={} 1={}", s, p);

    let e = Prob::from_prob(0.0);
    println!("0+0={} {}", e, e + e);

    println!("test passed");
}

pub fn test2() {
    let _x = bio::stats::LogProb::from(bio::stats::Prob(0.5));
    // XXX cannot multiply bio::stats::LogProb
    // println!("x*x={}", x * x);
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
    fn test_reflect() {
        let x = Prob::from_prob(1.0);
        let e = Prob::from_prob(0.0);
        assert_relative_eq!((x + e).0, x.0);
        assert_relative_eq!((x * e).0, e.0);
    }
    #[test]
    fn test_zero() {
        let zero = Prob::from_prob(0.0);
        println!("{:?}", zero);
        assert!(zero.is_zero());
        let nonzero = Prob::from_prob(0.00001);
        assert!(!nonzero.is_zero());
    }
}
