use crate::kmer::kmer::{tailing_kmers, Kmer};
use std::collections::HashMap;

pub trait DBG {
    fn new() -> Self;
    fn add(&mut self, kmer: Kmer, copy_num: u32);
    // fn update(self, kmer: Kmer, copy_num: u32);
    // fn remove(self, kmer: Kmer);
    fn find(&self, kmer: &Kmer) -> Option<u32>;
    fn is_exists(&self, kmer: &Kmer) -> bool;
    fn kmers(&self) -> Vec<Kmer>;
    // fn kmers_and_copy_nums(&self) -> (Vec<Kmer>, Vec<u32>);
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
}

pub struct DbgHash {
    store: HashMap<Kmer, u32>,
}

impl DBG for DbgHash {
    // implementation using hashmap as a store
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
    fn kmers(&self) -> Vec<Kmer> {
        self.store.keys().map(|x| x.clone()).collect()
    }
}

impl std::fmt::Display for DbgHash {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // generate dot graph file
        // digraph dbg {
        //   AATAT -> ATTTAT;
        // }
        writeln!(f, "digraph dbg {{");
        for kmer in self.kmers().iter() {
            // for node
            let copy_num = self.find(kmer).unwrap();
            writeln!(f, "\t{} [label=\"{} x{}\"];", kmer, kmer, copy_num);
            // for edges
            for child in self.childs(kmer).iter() {
                writeln!(f, "\t{} -> {};", kmer, child);
            }
        }
        writeln!(f, "}}");
        Ok(())
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
    let childs = d.childs(&Kmer::from(b"ATCG"));
    for child in childs {
        println!("child in store: {}", child);
    }
    let parents = d.parents(&Kmer::from(b"TCGT"));
    for parent in parents {
        println!("parent in store: {}", parent);
    }
    let parents = d.parents(&Kmer::from(b"ATCG"));
    for parent in parents {
        println!("parent in store: {}", parent);
    }
    for kmer in d.kmers().iter() {
        println!("kmer in store {}", kmer);
    }
    println!(":::DBG OUT:::\n{}", d);
}
