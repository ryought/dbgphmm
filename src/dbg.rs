use crate::hmm::base::Node;
use crate::kmer::kmer::{tailing_kmers, Kmer};
use crate::prob::Prob;
use arrayvec::ArrayVec;
use fnv::FnvHashMap as HashMap;
use log::{debug, info, warn};
use std::fmt::Write as FmtWrite;
// use ahash::AHashMap as HashMap;
// use std::collections::HashMap;

pub trait DBG {
    fn new() -> Self;
    fn add(&mut self, kmer: Kmer, copy_num: u32);
    // fn update(self, kmer: Kmer, copy_num: u32);
    // fn remove(self, kmer: Kmer);
    fn find(&self, kmer: &Kmer) -> u32;
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
    fn childs_with_copy_number(&self, kmer: &Kmer) -> Vec<(Kmer, u32)> {
        kmer.childs()
            .into_iter()
            .map(|child| {
                let copy_num = self.find(&child);
                (child, copy_num)
            })
            .filter(|(_, copy_num)| *copy_num > 0)
            .collect()
    }
    fn parents_with_copy_number(&self, kmer: &Kmer) -> Vec<(Kmer, u32)> {
        kmer.parents()
            .into_iter()
            .map(|parent| {
                let copy_num = self.find(&parent);
                (parent, copy_num)
            })
            .filter(|(_, copy_num)| *copy_num > 0)
            .collect()
    }
    fn childs_with_trans_prob(&self, kmer: &Kmer) -> Vec<(Kmer, Prob)> {
        let childs_with_cn = self.childs_with_copy_number(kmer);
        let sum_cn: u32 = childs_with_cn.iter().map(|(kmer, cn)| cn).sum();
        childs_with_cn
            .into_iter()
            .map(|(kmer, cn)| (kmer, Prob::from_prob(f64::from(cn) / f64::from(sum_cn))))
            .collect()
    }
    fn foremost_kmers(&self) -> Vec<Kmer> {
        // kmers with no parents
        self.kmers()
            .iter()
            .filter(|kmer| self.parents(kmer).len() == 0)
            .map(|kmer| kmer.clone())
            .collect()
    }
    fn vectorize(
        &self,
    ) -> (
        Vec<Kmer>,
        Vec<ArrayVec<Node, 4>>,
        Vec<ArrayVec<Node, 4>>,
        Vec<ArrayVec<Prob, 4>>,
    ) {
        // assign an index to each kmers
        // and returns (kmers, copy_nums, childs, parents, trans_probs)
        let kmers = self.kmers();
        let mut ids: HashMap<Kmer, usize> = HashMap::default();
        for (i, kmer) in kmers.iter().enumerate() {
            ids.insert(kmer.clone(), i);
        }

        let mut childs: Vec<ArrayVec<Node, 4>> = Vec::new();
        let mut trans_probs: Vec<ArrayVec<Prob, 4>> = Vec::new();
        let mut parents: Vec<ArrayVec<Node, 4>> = Vec::new();
        for kmer in kmers.iter() {
            let childs_with_tp = self.childs_with_trans_prob(kmer);
            childs.push(
                childs_with_tp
                    .iter()
                    .map(|(kmer, _)| Node(*ids.get(kmer).unwrap()))
                    .collect(),
            );
            trans_probs.push(childs_with_tp.iter().map(|(_, tp)| *tp).collect());
            parents.push(
                self.parents(kmer)
                    .iter()
                    .map(|kmer| Node(*ids.get(kmer).unwrap()))
                    .collect(),
            );
        }
        (kmers, childs, parents, trans_probs)
    }
    fn is_copy_number_consistent(&self) -> bool {
        // for all kmers
        for kmer in self.kmers().iter() {
            debug!("consistency checking {}", kmer);
            // check if sum_cn(childs) == sum_cn(siblings)
            // siblings = parents of a child
            let childs_with_cn = self.childs_with_copy_number(kmer);
            if childs_with_cn.len() > 0 {
                let sum_cn_childs: u32 = childs_with_cn.iter().map(|(_, cn)| cn).sum();
                let (child, _) = childs_with_cn.first().unwrap();
                let siblings_with_cn = self.parents_with_copy_number(child);
                let sum_cn_siblings: u32 = siblings_with_cn.iter().map(|(_, cn)| cn).sum();
                if sum_cn_childs != sum_cn_siblings {
                    return false;
                }
            }
        }
        true
    }
    fn augment_edge_kmers(&mut self) {
        // add `NATC`, `NNAT`, `NNNA` for foremost kmer `ATCG`
        for kmer in self.foremost_kmers().iter() {
            let cn = self.find(kmer);
            for new_kmer in tailing_kmers(kmer).into_iter() {
                self.add(new_kmer, cn);
            }
        }
    }
    fn as_dot(&self) -> String {
        // generate dot graph file
        // digraph dbg {
        //   AATAT -> ATTTAT;
        // }
        let mut s = String::new();
        writeln!(&mut s, "digraph dbg {{");
        writeln!(&mut s, "node [fontsize = 6, shape = box];");
        writeln!(&mut s, "edge [fontsize = 6];");
        for kmer in self.kmers().iter() {
            // for node
            let copy_num = self.find(kmer);
            writeln!(&mut s, "\t{} [label=\"{} x{}\"];", kmer, kmer, copy_num);
            // for edges
            for (child, p) in self.childs_with_trans_prob(kmer).iter() {
                writeln!(&mut s, "\t{} -> {} [label=\"{}\"];", kmer, child, p);
            }
        }
        writeln!(&mut s, "}}");
        s
    }
    fn as_degree_stats(&self) {
        let mut in_degs: [u32; 4] = [0; 4];
        let mut out_degs: [u32; 4] = [0; 4];
        for kmer in self.kmers().iter() {
            let in_deg = self.parents(kmer).len();
            let out_deg = self.childs(kmer).len();
            in_degs[in_deg] += 1;
            out_degs[out_deg] += 1;
        }
        println!("degree stats");
        for i in 0..4 {
            println!("in_deg={}: {}", i, in_degs[i]);
            println!("out_deg={}: {}", i, out_degs[i]);
        }
    }
}

pub struct DbgHash {
    store: HashMap<Kmer, u32>,
}

impl DBG for DbgHash {
    // implementation using hashmap as a store
    fn new() -> DbgHash {
        let h: HashMap<Kmer, u32> = HashMap::default();
        DbgHash { store: h }
    }
    fn add(&mut self, kmer: Kmer, copy_num: u32) {
        let copy_num_old = self.find(&kmer);
        self.store.insert(kmer, copy_num + copy_num_old);
    }
    fn find(&self, kmer: &Kmer) -> u32 {
        match self.store.get(&kmer) {
            Some(&copy_num) => copy_num,
            _ => 0,
        }
    }
    fn is_exists(&self, kmer: &Kmer) -> bool {
        self.store.contains_key(kmer)
    }
    fn kmers(&self) -> Vec<Kmer> {
        self.store.keys().map(|x| x.clone()).collect()
    }
}

impl DbgHash {
    pub fn from(kmers: Vec<Kmer>, copy_nums: Vec<u32>) -> DbgHash {
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
        Kmer::from(b"GGAC"),
        Kmer::from(b"TGAC"),
        Kmer::from(b"AGAC"),
        Kmer::from(b"GACT"),
        /*
        Kmer::from(b"TTCG"),
        Kmer::from(b"TCGT"),
        Kmer::from(b"TCGA"),
        */
    ];
    let copy_nums: Vec<u32> = vec![3, 2, 3, 3, 1];
    let mut d = DbgHash::from(kmers, copy_nums);
    d.augment_edge_kmers();
    println!("{}", d.as_dot());
}
