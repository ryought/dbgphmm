use criterion::{
    black_box, criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion,
    PlotConfiguration,
};
use std::{thread, time};
extern crate dbgphmm;
use dbgphmm::kmer::counter;

fn fibonacci(n: u64) -> u64 {
    match n {
        0 => 1,
        1 => 1,
        n => fibonacci(n - 1) + fibonacci(n - 2),
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    /*
    c.bench_function("count w", |b| {
        let s = counter::get_random_vec(1_000_000);
        b.iter(|| {
            // example 1
            // thread::sleep(time::Duration::from_secs(1));
            // fibonacci(black_box(20))
            //
            // example 2
            // let a = Kmer::from(black_box(b"ATCGATTAG"));
            // let b = Kmer::from(black_box(b"ATCGATTAG"));
            // assert_eq!(a, b);
            //
            // example 3
            let _h = counter::count(&s, 8);
        })
    });
    c.bench_function("count k20", |b| {
        b.iter(|| {
            let _h = counter::count(&s, 20);
        })
    });
    c.bench_function("count_my_kmer k20", |b| {});
    */

    // for k in [2, 4, 6, 8].iter() {
    //     group.bench_with_input(BenchmarkId::from_parameter(k), k, |b, &k| {
    //         b.iter(|| {
    //             let _h = counter::count(&s, k);
    //         });
    //     });
    // }

    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);

    let mut group = c.benchmark_group("count vs count_my_kmer");
    group.plot_config(plot_config);
    for len in [100, 10_000, 1_000_000].iter() {
        let s = counter::get_random_vec(*len);
        group.bench_with_input(BenchmarkId::new("count", len), len, |b, &len| {
            b.iter(|| {
                let _h = counter::count(&s, 8);
            });
        });
        group.bench_with_input(BenchmarkId::new("count_my_kmer", len), len, |b, &len| {
            b.iter(|| {
                let _h = counter::count_my_kmer(&s, 8);
            });
        });
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
