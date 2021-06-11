use dbgphmm::*;
use std::io::prelude::*;

fn test() {
    // generics
    let l1 = vec![34, 50, 25];
    println!("l1: {}", hoge::fuga::largest(&l1));

    // string
    hoge::mystring::test1();

    // DNA
    /*
    let x = hoge::seq::MyDNA::A;
    let c = hoge::seq::show(x);
    let c2 = hoge::seq::show(hoge::seq::complement(x));
    println!("c = {}, c2 = {}", c, c2);
    */

    test_struct::hoge::test1();
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    // let config = parse_config(&args);
    let config = kmer::counter::Config::new(&args);
    // run(config);
    // test5();
    // prob::test();
    // seq::test();
    // counter::run_counter(config);

    // my_vec::test();
    // seq::test2();
    // linkedlist::test2();

    // vec_of_vec::test();
    // kmer::test();

    // kmer::counter::test_counter();
    // sleeper::sleep();
    hmm::test();
    prob::test();
}

fn run2(config: kmer::counter::Config) {
    let mut f = std::fs::File::open(config.filename).expect("file not found");
    let mut contents = String::new();
    f.read_to_string(&mut contents).expect("cannot read file");
    println!("text: {}", contents);
}

fn length(seq: &[u8]) -> usize {
    // seq.len()
    let mut n = 0;
    for s in seq.iter() {
        n += 1;
    }
    n
}

use std::collections::HashMap;
fn kmer_count(seq: &[u8], k: usize) -> HashMap<&[u8], usize> {
    let mut count = HashMap::new();
    for i in 0..10 {
        count.insert(&seq[i..i + k], 1);
    }
    count
}

#[derive(Debug)]
enum DNA {
    A,
    C,
    G,
    T,
    N,
}
#[derive(Debug)]
struct Seq {
    array: Vec<DNA>,
}
/*
impl Seq {
    fn from_slice() -> Seq {
    }
}
*/
fn test5() {
    let s1 = Seq {
        array: vec![DNA::A, DNA::C, DNA::T],
    };
    println!("seq {:?}", s1)
}
