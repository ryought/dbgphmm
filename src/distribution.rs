use crate::prob::Prob;
use libm::{erf, sqrt};
use std::f64::consts::PI;

/// TODO should compute in log-space
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

///
/// calculate P(X=x) according to X ~ Normal(mu, sigma)
///
pub fn normal(x: f64, mu: f64, sigma: f64) -> Prob {
    let sigma_2 = sigma.powi(2);
    Prob::from_log_prob((-0.5 * (2.0 * PI * sigma_2).ln()) - ((x - mu).powi(2) / (2.0 * sigma_2)))
}

/// Convert base-coverage into kmer-coverage
///
/// Given base-coverage=C read-length=L error-rate=p, kmer-coverage is `C * (L-k)/L * (1-p)^k`
///
pub fn kmer_coverage(k: usize, read_length: usize, base_coverage: f64, p_error: Prob) -> f64 {
    let ratio_length = (read_length - k) as f64 / read_length as f64;
    let ratio_error = (1.0 - p_error.to_value()).powi(k as i32);
    base_coverage * ratio_length * ratio_error
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prob::p;

    #[test]
    fn normal_distribution() {
        // TODO WIP
        // check the distribution
        // check if no floating-point error happened
        for x in 0..200 {
            println!("{}", normal_bin(x, 100, 10));
        }

        // check if the sum is 1
        let s: Prob = (0..200).map(|x| normal_bin(x, 100, 10)).sum();
        println!("{}", s);
    }

    #[test]
    fn normal_2() {
        let mut z = Prob::zero();
        for x in -10..10 {
            let p = normal(x as f64, 0.0, 1.0);
            println!("p({})={}", x, p);
            z += p;
        }
        assert!((normal(0.0, 0.0, 1.0).to_value() - 0.3989).abs() < 0.0001);
        assert!((normal(1.0, 0.0, 1.0).to_value() - 0.2420).abs() < 0.0001);
        println!("z={}", z);
        assert!((z.to_value() - 1.0).abs() < 0.000001);
    }

    #[test]
    fn kmer_coverage_test() {
        println!("{}", kmer_coverage(16, 100, 10.0, p(0.0)));
        println!("{}", kmer_coverage(16, 10_000, 10.0, p(0.0)));
        println!("{}", kmer_coverage(40, 10_000, 10.0, p(0.00_1))); // 0.1% 10k read
        println!("{}", kmer_coverage(32, 50, 20.0, p(0.00_1))); // small example
    }
}
