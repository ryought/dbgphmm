use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::{thread, time};

fn fibonacci(n: u64) -> u64 {
    match n {
        0 => 1,
        1 => 1,
        n => fibonacci(n - 1) + fibonacci(n - 2),
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("fib 20", |b| {
        b.iter(|| {
            thread::sleep(time::Duration::from_secs(1));
            fibonacci(black_box(20))
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
