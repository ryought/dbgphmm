///
/// probability calculation
/// implements logaddexp
///

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Prob {
    log_value: f64,
}

impl Prob {
    fn from_prob(value: f64) -> Prob {
        Prob {
            log_value: value.ln(),
        }
    }
    fn from_log_prob(log_value: f64) -> Prob {
        Prob { log_value }
    }
    fn to_value(self) -> f64 {
        self.log_value.exp()
    }
    fn to_log_value(self) -> f64 {
        self.log_value
    }
}

// display
impl std::fmt::Display for Prob {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.log_value)
    }
}

// operators
impl std::ops::Add for Prob {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        let x = self.log_value;
        let y = other.log_value;
        if x > y {
            Prob {
                log_value: x + ((y - x).exp() + 1.0).ln(),
            }
        } else {
            Prob {
                log_value: y + ((x - y).exp() + 1.0).ln(),
            }
        }
    }
}
impl std::ops::Mul for Prob {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Prob {
            log_value: self.log_value + other.log_value,
        }
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

pub fn test() {
    let z0 = Prob::from_prob(0.0);
    let z1 = Prob::from_prob(1.0);
    let x = Prob::from_prob(0.3);
    let y = Prob::from_prob(0.3);
    let A = Prob::from_prob(0.6);
    let a = x + y;
    let B = Prob::from_prob(0.09);
    let b = x * y;
    println!(
        "{} {} {} {} {} {} {}",
        x.to_value(),
        y.to_value(),
        a.to_value(),
        b.to_value(),
        z0.log_value,
        z1.log_value,
        (z0 + z1).to_value(),
    );
    assert_abs_diff_eq!(B.log_value, b.log_value);

    let x = Prob::from_prob(0.3);
    let e = Prob::from_prob(0.0);
    assert_relative_eq!((x + e).log_value, x.log_value);
    assert_relative_eq!((x * e).log_value, e.log_value);
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
    assert_relative_eq!(x.log_value, y.log_value);
    println!(
        "{} {} {} {}",
        x.to_value(),
        s.to_value(),
        y.to_value(),
        p.to_value()
    );

    println!("test passed");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_id() {
        let x = Prob::from_prob(0.3);
        let e = Prob::from_prob(0.0);
        assert_relative_eq!((x + e).log_value, x.log_value);
        assert_relative_eq!((x * e).log_value, e.log_value);
    }
    #[test]
    fn test_id_id() {
        let x = Prob::from_prob(1.0);
        let e = Prob::from_prob(0.0);
        assert_relative_eq!((x + e).log_value, x.log_value);
        assert_relative_eq!((x * e).log_value, e.log_value);
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
}
