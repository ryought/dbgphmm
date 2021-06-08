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

pub fn test() {
    let x = Prob::from_prob(0.003);
    let y = Prob::from_prob(0.003);
    let a = x * y;
    let b = x + y;
    println!(
        "{} {} {} {}",
        x.to_value(),
        y.to_value(),
        a.to_value(),
        b.to_value()
    );
    // println!("{}", x == y);
}
