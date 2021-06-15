struct Kmer {
    hoge: u8,
}

impl Kmer {
    fn new() -> Self {
        Self { hoge: 10 }
    }
}

fn my_iter(kmers: Vec<Kmer>) -> impl std::iter::Iterator {
    kmers.iter()
}
