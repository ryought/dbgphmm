#[cfg(test)]
mod tests {
    //
    // speed benchmarks
    //
    use test::Bencher;

    //
    // common settings
    //
    use crate::genome;
    use crate::prelude::*;
    use crate::utils::timer;
    fn g1m() -> Genome {
        // unit 1kb x 1000 times = 1MB diploid
        let genome = genome::tandem_repeat_polyploid_with_unique_homo_ends(
            1_000, 1_000, 0, 0.0, 0, 1_000, 2, 0.01, 0,
        );
        genome
    }

    #[bench]
    #[ignore = "~7s in release"]
    fn bench_g1m_k100_construct(b: &mut Bencher) {
        let genome = g1m();
        genome.show();

        // ~200ms
        let k = 100;
        let (dbg, time): (SimpleDbg<VecKmer>, _) =
            timer(|| SimpleDbg::from_styled_seqs(k, &genome));
        println!("time={}", time);

        b.iter(|| {
            // let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_styled_seqs(k, &genome);
            let k = 100;
            let (dbg, time): (SimpleDbg<VecKmer>, _) =
                timer(|| SimpleDbg::from_styled_seqs(k, &genome));
            println!("time={}", time);
            dbg
        });
    }

    #[bench]
    #[ignore = "~3s in release"]
    fn bench_g1m_k100_clone(b: &mut Bencher) {
        let genome = g1m();
        let k = 100;
        let (dbg, time): (SimpleDbg<VecKmer>, _) =
            timer(|| SimpleDbg::from_styled_seqs(k, &genome));
        println!("time={}", time);
        b.iter(|| dbg.clone())
    }
}
