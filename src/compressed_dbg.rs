use crate::cycles;
use crate::cycles::CycleDirection;
use crate::dbg::{DbgHash, DBG};
use crate::distribution::normal_bin;
use crate::graph::Node;
use crate::kmer::kmer::{linear_seq_to_kmers, null_kmer, Kmer};
use crate::prob::Prob;
use fnv::FnvHashMap as HashMap;
use histo::Histogram;
use log::{debug, info, warn};
use std::fmt::Write as FmtWrite;

/// indexed and compressed dbg
/// O(1) access to childs, parents, emissions of index
/// back to kmer content if necessary
/// and it stores indexed cycles
pub struct CompressedDBG {
    k: usize,
    n_kmers: usize,
    childs: Vec<Vec<Node>>,
    parents: Vec<Vec<Node>>,
    kmers: Vec<Kmer>,
    ids: HashMap<Kmer, Node>,
    emissions: Vec<u8>,
    cycles: Vec<Vec<(Node, CycleDirection)>>,
}

impl CompressedDBG {
    /// Construct CompressedDBG from DbgHash
    /// - assign an index to each kmers and convert hashmap into vector
    /// - calc cycles and store in vector
    pub fn from(dbg: &DbgHash, k: usize) -> CompressedDBG {
        // assign an index to each kmers
        let kmers: Vec<Kmer> = dbg.kmers();
        let n_kmers = kmers.len();
        let mut ids: HashMap<Kmer, Node> = HashMap::default();
        for (i, kmer) in kmers.iter().enumerate() {
            ids.insert(kmer.clone(), Node(i));
        }

        let mut childs: Vec<Vec<Node>> = Vec::new();
        let mut parents: Vec<Vec<Node>> = Vec::new();
        let mut emissions: Vec<u8> = Vec::new();
        for kmer in kmers.iter() {
            childs.push(
                dbg.childs(kmer)
                    .iter()
                    .map(|kmer| *ids.get(kmer).unwrap())
                    .collect(),
            );
            parents.push(
                dbg.parents(kmer)
                    .iter()
                    .map(|kmer| *ids.get(kmer).unwrap())
                    .collect(),
            );
            emissions.push(kmer.last());
        }

        let root = null_kmer(k - 1);
        let s = cycles::DbgTree::new(dbg, &root);

        let cycles: Vec<Vec<(Node, CycleDirection)>> = s
            .cycle_keys()
            .iter()
            .map(|key| {
                s.cycle_components(key)
                    .iter()
                    .map(|(kmer, dir)| (*ids.get(kmer).unwrap(), *dir))
                    .collect()
            })
            .collect();

        CompressedDBG {
            k,
            n_kmers,
            childs,
            parents,
            kmers,
            ids,
            emissions,
            cycles,
        }
    }
    // TODO add construction method from k-1mer compresseddbg
    pub fn n_kmers(&self) -> usize {
        self.n_kmers
    }
    pub fn iter_nodes(&self) -> impl std::iter::Iterator<Item = Node> {
        (0..self.n_kmers).map(|i| Node(i))
    }
    pub fn childs(&self, v: &Node) -> &[Node] {
        &self.childs[v.0]
    }
    pub fn parents(&self, v: &Node) -> &[Node] {
        &self.parents[v.0]
    }
    pub fn is_adjacent(&self, v: &Node, w: &Node) -> bool {
        self.childs(v).iter().any(|u| *u == *w)
    }
    pub fn emission(&self, v: &Node) -> u8 {
        self.emissions[v.0]
    }
    pub fn is_emitable(&self, v: &Node) -> bool {
        self.emission(v) != b'N'
    }
    pub fn kmer(&self, v: &Node) -> &Kmer {
        &self.kmers[v.0]
    }
    pub fn n_cycles(&self) -> usize {
        self.cycles.len()
    }
    pub fn cycle_components(&self, cycle_id: usize) -> &[(Node, CycleDirection)] {
        &self.cycles[cycle_id]
    }
    /// Get a vec of NNNNN{ACGT} in cdbg
    pub fn heads(&self) -> Vec<Node> {
        null_kmer(self.k - 1)
            .succs()
            .iter()
            .filter_map(|kmer| self.ids.get(&kmer))
            .map(|&v| v)
            .collect()
    }
    /// Calc transition probabilities from copy numbers of each kmers
    /// given
    /// childs[i] = [0, 1, 2, 3]
    /// this computes
    /// trans_probs[i] = [p(i->0), p(i->1), p(i->2), p(i->3)]
    pub fn copy_num_to_trans_prob(&self, copy_nums: &[u32]) -> Vec<Vec<Prob>> {
        assert_eq!(self.n_kmers(), copy_nums.len());
        self.iter_nodes()
            .map(|v| {
                let cns: Vec<u32> = self
                    .childs(&v)
                    .iter()
                    .map(|&w| {
                        if self.is_emitable(&w) {
                            copy_nums[w.0]
                        } else {
                            0
                        }
                    })
                    .collect();
                let total_cn: u32 = cns.iter().sum();
                cns.iter()
                    .map(|&cn| {
                        if self.is_emitable(&v) && total_cn > 0 {
                            Prob::from_prob(f64::from(cn) / f64::from(total_cn))
                        } else {
                            Prob::from_prob(0.0)
                        }
                    })
                    .collect()
            })
            .collect()
    }
    /// check copy_nums consistency by
    /// (total cn of parents) == (total cn of childs of one of the parents)
    pub fn is_consistent_copy_num(&self, copy_nums: &[u32]) -> bool {
        if self.n_kmers() != copy_nums.len() {
            warn!("copy_nums shape mismatch");
            return false;
        }
        self.iter_nodes().all(|v| {
            let parents = self.parents(&v);
            if parents.len() == 0 {
                false
            } else {
                let siblings = self.childs(&parents[0]);
                let total_parents: u32 = parents.iter().map(|v| copy_nums[v.0]).sum();
                let total_siblings: u32 = siblings.iter().map(|v| copy_nums[v.0]).sum();
                total_parents == total_siblings
            }
        })
    }
    pub fn total_emitable_copy_num(&self, copy_nums: &[u32]) -> u32 {
        (0..self.n_kmers())
            .map(|i| Node(i))
            .filter(|&v| self.is_emitable(&v))
            .map(|v| copy_nums[v.0])
            .sum()
    }
    /// cycle is a sequence of (node, direction)
    /// when is_up=true, node with direction=Forward will get +1.
    /// all nodes should be >=0 after this update
    pub fn is_acceptable(&self, copy_nums: &[u32], cycle_id: usize, is_up: bool) -> bool {
        self.cycle_components(cycle_id).iter().all(|(v, dir)| {
            if is_up && *dir == CycleDirection::Reverse {
                copy_nums[v.0] > 0
            } else if !is_up && *dir == CycleDirection::Forward {
                copy_nums[v.0] > 0
            } else {
                true
            }
        })
    }
    pub fn update_by_cycle(&self, copy_nums: &[u32], cycle_id: usize, is_up: bool) -> Vec<u32> {
        let mut new_copy_nums = copy_nums.to_vec();
        for (v, dir) in self.cycle_components(cycle_id).iter() {
            if is_up {
                match dir {
                    CycleDirection::Forward => new_copy_nums[v.0] += 1,
                    CycleDirection::Reverse => new_copy_nums[v.0] -= 1,
                }
            } else {
                match dir {
                    CycleDirection::Forward => new_copy_nums[v.0] -= 1,
                    CycleDirection::Reverse => new_copy_nums[v.0] += 1,
                }
            }
        }
        new_copy_nums
    }
    pub fn is_all_cycle_consistent(&self, init_copy_nums: &[u32]) {
        for i in 0..self.n_cycles() {
            let x = self.is_acceptable(init_copy_nums, i, true);
            let y = self.is_acceptable(init_copy_nums, i, false);
            println!("i={}, x={}, y={}", i, x, y);
        }
    }
    pub fn from_seqs(seqs: Vec<Vec<u8>>, k: usize) -> (CompressedDBG, Vec<u32>) {
        let dbg = DbgHash::from_seqs(seqs, k);
        let cdbg = CompressedDBG::from(&dbg, k);
        let copy_nums: Vec<u32> = cdbg
            .iter_nodes()
            .map(|v| {
                let kmer = cdbg.kmer(&v);
                dbg.find(kmer)
            })
            .collect();
        (cdbg, copy_nums)
    }
    /// calculate true copy nums from seq
    /// Main use case is to determine the true copy_nums on cdbg constructed from reads
    /// by using reference information
    pub fn true_copy_nums_from_seqs(&self, seqs: Vec<Vec<u8>>, k: usize) -> Option<Vec<u32>> {
        let mut copy_nums: Vec<u32> = vec![0; self.n_kmers()];
        for seq in seqs.iter() {
            for kmer in linear_seq_to_kmers(seq, k) {
                match self.ids.get(&kmer) {
                    Some(&v) => copy_nums[v.0] += 1,
                    None => {
                        warn!("true kmer {} not found in cdbg", kmer);
                        return None;
                    }
                }
            }
        }
        Some(copy_nums)
    }
    pub fn check_kmer_existence(&self, seqs: Vec<Vec<u8>>, k: usize) {
        let mut t: u32 = 0;
        let mut f: u32 = 0;
        for seq in seqs.iter() {
            for kmer in linear_seq_to_kmers(seq, k) {
                match self.ids.get(&kmer) {
                    Some(&v) => t += 1,
                    None => f += 1,
                }
            }
        }
        warn!("exist={} not_exist={}", t, f);
    }
    /// prior score of this
    /// Assuming genome size ~ Normal(ave, std)
    pub fn prior_score(&self, copy_nums: &[u32], ave_size: u32, std_size: u32) -> Prob {
        normal_bin(self.total_emitable_copy_num(copy_nums), ave_size, std_size)
    }
    /// Graphviz dot format
    pub fn as_dot(&self) -> String {
        let mut s = String::new();
        writeln!(&mut s, "digraph cdbg {{");
        for v in self.iter_nodes() {
            // for node
            let kmer = self.kmer(&v);
            writeln!(&mut s, "\t{} [label=\"{}\"];", v.0, kmer);
            // for edges
            for w in self.childs(&v).iter() {
                // writeln!(&mut s, "\t{} -> {} [label=\"{}\"];", v.0, w.0);
                writeln!(&mut s, "\t{} -> {};", v.0, w.0);
            }
        }
        writeln!(&mut s, "}}");
        s
    }
    /// Graphviz dot format
    pub fn as_dot_with_cycle(&self, cycle_id: usize) -> String {
        let cycle = self.cycle_components(cycle_id);
        let mut s = String::new();
        writeln!(&mut s, "digraph cdbg {{");
        for v in self.iter_nodes() {
            // for node
            let kmer = self.kmer(&v);
            if cycle.iter().any(|&(w, _)| w == v) {
                writeln!(&mut s, "\t{} [label=\"{}\" color=red];", v.0, kmer);
            } else {
                writeln!(&mut s, "\t{} [label=\"{}\"];", v.0, kmer);
            }
            // for edges
            for w in self.childs(&v).iter() {
                // writeln!(&mut s, "\t{} -> {} [label=\"{}\"];", v.0, w.0);
                writeln!(&mut s, "\t{} -> {};", v.0, w.0);
            }
        }
        writeln!(&mut s, "}}");
        s
    }
    /// dot with copy_number info on nodes
    pub fn as_dot_with_copy_nums(&self, copy_nums: &[u32]) -> String {
        let mut s = String::new();
        writeln!(&mut s, "digraph cdbg {{");
        for v in self.iter_nodes() {
            // for node
            let kmer = self.kmer(&v);
            let copy_num = copy_nums[v.0];
            writeln!(&mut s, "\t{} [label=\"{} x{}\"];", v.0, kmer, copy_num);
            // for edges
            for w in self.childs(&v).iter() {
                writeln!(&mut s, "\t{} -> {};", v.0, w.0);
            }
        }
        writeln!(&mut s, "}}");
        s
    }

    pub fn as_cycle_stats(&self) -> String {
        let mut s = String::new();
        writeln!(&mut s, "#cycles={}", self.n_cycles());
        let mut histogram = Histogram::with_buckets(10);
        for i in 0..self.n_cycles() {
            histogram.add(self.cycle_components(i).len() as u64);
        }
        writeln!(&mut s, "{}", histogram);
        s
    }

    pub fn as_degree_stats(&self) -> String {
        let mut s = String::new();

        let mut in_degs: [u32; 6] = [0; 6];
        let mut out_degs: [u32; 6] = [0; 6];
        for v in self.iter_nodes() {
            let in_deg = self.parents(&v).len();
            let out_deg = self.childs(&v).len();
            in_degs[in_deg] += 1;
            out_degs[out_deg] += 1;
        }
        writeln!(&mut s, "#kmer: {}", self.n_kmers());
        writeln!(&mut s, "indeg");
        for i in 0..6 {
            writeln!(&mut s, "{}:\t{}", i, in_degs[i]);
        }
        writeln!(&mut s, "outdeg");
        for i in 0..6 {
            writeln!(&mut s, "{}:\t{}", i, out_degs[i]);
        }
        s
    }
}

pub fn test() {
    let mut d = DbgHash::new();
    // ATCGATTCGAT
    d.add_seq(b"ATCGATTCGATCGATTCGATAGATCG", 8);
    // d.add_seq(b"ATCTCGATCTTGATAGATCG", 8);
    // println!("{}", d.as_dot());
    let cdbg = CompressedDBG::from(&d, 8);
    let copy_nums_true: Vec<u32> = cdbg
        .iter_nodes()
        .map(|v| {
            let kmer = cdbg.kmer(&v);
            d.find(kmer)
        })
        .collect();
    println!("copy-nums {:?}", copy_nums_true);

    let p = cdbg.is_consistent_copy_num(&copy_nums_true);
    println!("copy-nums-consistency {:?}", p);

    println!("#cycles: {}", cdbg.n_cycles());
    for i in 0..cdbg.n_cycles() {
        println!("cycle#{}: {:?}", i, cdbg.cycle_components(i));
    }

    /*
    for (i, v) in cdbg.iter_nodes().enumerate() {
        println!(
            "node #{}: {:?} childs:{:?} parents:{:?}",
            i,
            v,
            cdbg.childs(v),
            cdbg.parents(v)
        );
    }
    println!("#cycles: {}", cdbg.n_cycles());
    for i in 0..cdbg.n_cycles() {
        println!("cycle#{}: {:?}", i, cdbg.cycle_components(i));
    }
    */

    // println!("{}", cdbg.as_dot_with_cycle(1));
}
