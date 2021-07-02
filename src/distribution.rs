use crate::prob::Prob;
use libm::{erf, sqrt};

pub fn normal_cdf(x: u32, mu: u32, sigma: u32) -> Prob {
    let f_x = x as f64;
    let f_mu = mu as f64;
    let f_sigma = sigma as f64;
    Prob::from_prob(0.5 * (1.0 + erf((f_x - f_mu) / sqrt(2.0 * f_sigma * f_sigma))))
}

pub fn normal_bin(x: u32, mu: u32, sigma: u32) -> Prob {
    let p0 = normal_cdf(x + 1, mu, sigma);
    let p1 = normal_cdf(x, mu, sigma);
    Prob::from_prob(p0.to_value() - p1.to_value())
}

pub fn test() {
    for x in 0..100 {
        println!("{}", normal_bin(x, 50, 10));
    }
}
