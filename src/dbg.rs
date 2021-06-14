use crate::kmer::kmer::{tailing_kmers, Kmer};
use std::collections::HashMap;

pub trait DBG {
    fn new() -> Self;
    fn add(&mut self, kmer: Kmer, copy_num: u32);
    // fn update(self, kmer: Kmer, copy_num: u32);
    // fn remove(self, kmer: Kmer);
    fn find(&self, kmer: &Kmer) -> Option<u32>;
    fn is_exists(&self, kmer: &Kmer) -> bool;
    fn childs(&self, kmer: &Kmer) -> Vec<Kmer>;
    fn parents(&self, kmer: &Kmer) -> Vec<Kmer>;
    fn as_dot(&self) -> String;
}

pub struct DbgHash {
    store: HashMap<Kmer, u32>,
}

impl DBG for DbgHash {
    // store functions
    fn new() -> DbgHash {
        let h: HashMap<Kmer, u32> = HashMap::new();
        DbgHash { store: h }
    }
    fn add(&mut self, kmer: Kmer, copy_num: u32) {
        self.store.insert(kmer, copy_num);
    }
    fn find(&self, kmer: &Kmer) -> Option<u32> {
        match self.store.get(&kmer) {
            Some(&copy_num) => Some(copy_num),
            _ => None,
        }
    }
    fn is_exists(&self, kmer: &Kmer) -> bool {
        self.store.contains_key(kmer)
    }
    fn childs(&self, kmer: &Kmer) -> Vec<Kmer> {
        kmer.childs()
            .into_iter()
            .filter(|child| self.is_exists(child))
            .collect()
    }
    fn parents(&self, kmer: &Kmer) -> Vec<Kmer> {
        kmer.parents()
            .into_iter()
            .filter(|parent| self.is_exists(parent))
            .collect()
    }
    // output related
    fn as_dot(&self) -> String {
        // generate dot graph file
        "hogehoge".to_string()
    }
}

impl DbgHash {
    fn from(kmers: Vec<Kmer>, copy_nums: Vec<u32>) -> DbgHash {
        let mut d = DbgHash::new();
        for (kmer, copy_num) in kmers.into_iter().zip(copy_nums.into_iter()) {
            d.add(kmer, copy_num);
        }
        d
    }
}

pub fn test() {
    let kmers: Vec<Kmer> = vec![
        Kmer::from(b"ATCG"),
        Kmer::from(b"TTCG"),
        Kmer::from(b"TCGT"),
        Kmer::from(b"TCGA"),
    ];
    let copy_nums: Vec<u32> = vec![1, 1, 1, 1];
    let d = DbgHash::from(kmers, copy_nums);
    println!("{}", d.as_dot());
    let childs = Kmer::from(b"ATCG").childs();
    for child in childs {
        println!("all child: {}", child);
    }
    let parents = Kmer::from(b"ATCG").parents();
    for parent in parents {
        println!("all parent: {}", parent);
    }
    let childs = d.childs(&Kmer::from(b"ATCG"));
    for child in childs {
        println!("child in store: {}", child);
    }
}
