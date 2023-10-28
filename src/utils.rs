use indicatif::ProgressStyle;

pub fn is_equal_as_set<T: PartialEq>(xs: &[T], ys: &[T]) -> bool {
    xs.iter().all(|x| ys.contains(x)) && ys.iter().all(|y| xs.contains(y))
}

///
/// Check if all items from iterator is the same value
/// If so, return the unique value.
/// Otherwise (including the iterator has no value) return none.
///
pub fn all_same_value<T: PartialEq + Clone, I: Iterator<Item = T>>(mut iter: I) -> Option<T> {
    match iter.next() {
        Some(first) => {
            if iter.all(|item| item == first) {
                Some(first)
            } else {
                None
            }
        }
        None => None,
    }
}

///
/// TODO
/// * change return type from Vec<T> into Option<Vec<T>>?
///
pub fn unwrap_all<T>(options: Vec<Option<T>>) -> Vec<T> {
    options
        .into_iter()
        .map(|option| option.expect(""))
        .collect()
}

///
/// Given a list `xs` of (Value: T, Key: K) get breakpoints that is key of xs[i-1] and xs[i] is different.
///
pub fn breakpoints<T, K: PartialEq>(xs: &[(T, K)]) -> Vec<usize> {
    let n = xs.len();
    let mut ret = Vec::new();
    for i in 0..n {
        let im1 = if i == 0 { n - 1 } else { i - 1 };
        if xs[im1].1 != xs[i].1 {
            ret.push(i);
        }
    }
    ret
}

use std::time::Instant;
///
/// measure time in milli-seconds (ms) of closure.
///
pub fn timer<F, T>(f: F) -> (T, u128)
where
    F: FnOnce() -> T,
{
    let start = Instant::now();
    let ret = f();
    let duration = start.elapsed();
    (ret, duration.as_millis())
}

///
/// measure time in micro seconds (us) of closure.
///
pub fn timer_us<F, T>(f: F) -> (T, u128)
where
    F: FnOnce() -> T,
{
    let start = Instant::now();
    let ret = f();
    let duration = start.elapsed();
    (ret, duration.as_micros())
}

///
/// get strings with repeated n-times space (' ').
///
pub fn spaces(n: usize) -> String {
    // old rust
    // std::iter::repeat(" ").take(n).collect::<String>()
    // new rust 1.16
    " ".repeat(n)
}

pub fn check_memory_usage() {
    jemalloc_ctl::epoch::advance().unwrap();
    let allocated = jemalloc_ctl::stats::allocated::read().unwrap();
    let resident = jemalloc_ctl::stats::resident::read().unwrap();
    eprintln!("[memory] {} / {}", allocated, resident);
}

///
/// get dbgphmm/resources directory
///
pub fn resource_dir() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("resources")
}

///
/// calculate `log(n!)` by `log(n) + log(n-1) + ... + log(1)`
///
pub fn log_factorial(n: usize) -> f64 {
    let mut ret = 0.0;
    for i in (1..=n).rev() {
        ret += (i as f64).ln();
    }
    ret
}

pub fn progress_common_style() -> ProgressStyle {
    ProgressStyle::with_template(
        "[{elapsed_precise}/{eta_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
    )
    .unwrap()
    .progress_chars("##-")
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_equal_as_set_test_1() {
        let xs = vec![0, 1, 2, 5, 3];
        let ys = vec![1, 5, 0, 3, 2];
        assert!(is_equal_as_set(&xs, &ys));

        let zs = vec![1, 5, 3, 2];
        assert!(!is_equal_as_set(&xs, &zs));
    }

    #[test]
    fn all_same_value_test() {
        {
            let xs = vec![0, 1, 2, 5, 3];
            assert_eq!(all_same_value(xs.iter()), None);
        }
        {
            let xs = vec![0, 0, 0];
            assert_eq!(all_same_value(xs.iter()), Some(&0));
        }
        {
            let xs = vec![0, 0, 1];
            assert_eq!(all_same_value(xs.iter()), None);
        }
        {
            let xs = vec![5];
            assert_eq!(all_same_value(xs.iter()), Some(&5));
        }
        {
            let xs: Vec<usize> = Vec::new();
            assert_eq!(all_same_value(xs.iter()), None);
        }
    }

    #[test]
    fn log_factorial_test() {
        println!("log(0!)={}", log_factorial(0));
        println!("log(1!)={}", log_factorial(1));
        println!("log(2!)={}", log_factorial(2));
        println!("log(4!)={}", log_factorial(4));

        for n in 0..10 {
            let t: usize = (1..=n).product();
            println!("n={} n!={}", n, t);
            let lf = log_factorial(n);
            let t_approx = lf.exp();
            println!("n!~{}", t_approx);
            let diff = (t_approx - t as f64).abs();
            println!("diff={}", diff);
            assert!(diff < 0.00001);
        }
    }
}
