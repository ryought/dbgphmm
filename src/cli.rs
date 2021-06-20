use crate::hmm::params::PHMMParams;
use crate::random_seq;

pub fn generate(length: usize, seed: u64) {
    let v = random_seq::generate(length, seed);
    println!(">randseq");
    println!("{}", std::str::from_utf8(&v).unwrap());
}

pub fn sample(dbg_fa: String, length: u32, n_reads: u32, param: PHMMParams) {
    println!("sampling {} {} {}", dbg_fa, length, n_reads);
    println!("with param {}", param);
}

/*
fn calc_prob() {
    let (kmers, copy_nums) = io::fasta::parse_kmers_and_copy_nums(&opts.dbg_fa, opts.k);
    info!("from dbg_fa #kmer:{}", kmers.len());
    let d = hmm::dbg::DbgPHMM::new(kmers, copy_nums).unwrap();
    let reads = io::fasta::parse_reads(&opts.reads_fa);
    let param = hmm::params::PHMMParams::default();
    let p = d.forward_prob(&param, &reads[0]);
    println!("forward prob : {}", p);
    // let p = d.backward(&param, &reads[0]);
    // println!("backward prob : {}", p[0].FM);
}
*/
