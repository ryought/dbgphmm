struct Kmer {
    id: u64,
    seq: Vec<u8>,
    count: u64,
}

impl Kmer {
    fn new() {}
    fn is_parent() {}
    fn is_child() {}
}

/// check flow
struct DeBruijnGraph {
    kmers: Vec<Kmer>,
}

impl DeBruijnGraph {
    fn parents() {}
    fn childs() {}
}
