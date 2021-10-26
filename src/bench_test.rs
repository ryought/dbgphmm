fn fibonacci(n: u64) -> u64 {
    match n {
        0 => 1,
        1 => 1,
        n => fibonacci(n - 1) + fibonacci(n - 2),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use test::Bencher;

    #[bench]
    fn bench_fibonacci(b: &mut Bencher) {
        b.iter(|| {
            let n = test::black_box(20);
            fibonacci(n)
        })
        /*
        b.iter(|| {
            let n = test::black_box(1000);
            (0..n).fold(0, |a, b| a ^ b)
        })
        */
    }
}
