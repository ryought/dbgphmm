use crate::cycles;
use crate::cycles::CycleDirection;
use crate::dbg::{DbgHash, DBG};
use crate::distribution::normal_bin;
use crate::graph::{Edge, Edgei, Edgep, IndexedDiGraph, Node, Pos};
use crate::io::cytoscape::Element;
use crate::kmer::kmer::{linear_seq_to_kmers, null_kmer, Kmer};
use crate::prob::Prob;
use crate::stats;
use fnv::FnvHashMap as HashMap;
use fnv::FnvHashSet as HashSet;
use histo::Histogram;
use itertools::iproduct;
use log::{debug, info, warn};
use std::fmt::Write as FmtWrite;

/// indexed and compressed dbg
/// O(1) access to childs, parents, emissions of index
/// back to kmer content if necessary
/// and it stores indexed cycles
#[derive(PartialEq, Debug)]
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
    pub fn k(&self) -> usize {
        self.k
    }
    pub fn is_valid_node(&self, v: &Node) -> bool {
        0 <= v.0 && v.0 < self.n_kmers
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
    pub fn id(&self, kmer: &Kmer) -> Option<Node> {
        match self.ids.get(kmer) {
            Some(&v) => Some(v),
            None => None,
        }
    }
    pub fn edgep_to_edgei(&self, e: &Edgep) -> Edgei {
        let v = e.source;
        let w = e.target;
        let i = self.childs(&v).iter().position(|u| u == &w).unwrap();
        Edgei::new(v, i)
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
    /// Get a vec of N{ACGT}^{k-1} in cdbg
    pub fn starting_kmers(&self) -> Vec<Node> {
        self.iter_nodes()
            .filter_map(|v| {
                let kmer = self.kmer(&v);
                if kmer.is_starting() {
                    Some(v)
                } else {
                    None
                }
            })
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
                // TODO check my own copy numbers > 0 ?
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
    pub fn update_cycle_vec_by_cycle(
        &self,
        cycle_vec: &[u32],
        cycle_id: usize,
        is_up: bool,
    ) -> Vec<u32> {
        let mut new_cycle_vec = cycle_vec.to_vec();
        if is_up {
            new_cycle_vec[cycle_id] += 1;
        } else {
            new_cycle_vec[cycle_id] -= 1;
        }
        new_cycle_vec
    }
    /// Convert copy_nums into cycle_vec
    /// For each cycle basis,
    /// (the amount of cycle basis in the cycle vec) == (copy_num of cycle key edge)
    /// and the total amount will be the cycle vec
    pub fn cycle_vec_from_copy_nums(&self, copy_nums: &[u32]) -> Vec<u32> {
        let mut cycle_vec: Vec<i32> = vec![0; self.n_cycles()];
        for cycle_id in 0..self.n_cycles() {
            let (node, _) = self.cycle_components(cycle_id).first().unwrap();
            let count = copy_nums[node.0];
            cycle_vec[cycle_id] += count as i32;
        }
        cycle_vec
            .iter()
            .map(|&x| {
                assert!(x >= 0);
                x as u32
            })
            .collect()
    }
    /// Convert cycle_vec into copy_nums
    /// Inverse function of cycle_vec_from_copy_nums
    pub fn copy_nums_from_cycle_vec(&self, cycle_vec: &[u32]) -> Vec<u32> {
        let mut copy_nums: Vec<i32> = vec![0; self.n_kmers()];
        for cycle_id in 0..self.n_cycles() {
            let count = cycle_vec[cycle_id];
            for (node, dir) in self.cycle_components(cycle_id).iter() {
                match dir {
                    CycleDirection::Forward => copy_nums[node.0] += count as i32,
                    CycleDirection::Reverse => copy_nums[node.0] -= count as i32,
                }
            }
        }
        copy_nums
            .iter()
            .map(|&x| {
                assert!(x >= 0);
                x as u32
            })
            .collect()
    }
    pub fn is_all_cycle_consistent(&self, init_copy_nums: &[u32]) {
        for i in 0..self.n_cycles() {
            let x = self.is_acceptable(init_copy_nums, i, true);
            let y = self.is_acceptable(init_copy_nums, i, false);
            println!("i={}, x={}, y={}", i, x, y);
        }
    }
    pub fn validate_cycles(&self, copy_nums_true: &[u32]) {
        for i in 0..self.n_cycles() {
            let is_a = self.is_acceptable(&copy_nums_true, i, true);
            let is_b = self.is_acceptable(&copy_nums_true, i, false);
            assert!(is_a);
            assert!(is_b);
            let new_a = self.update_by_cycle(&copy_nums_true, i, true);
            let new_b = self.update_by_cycle(&copy_nums_true, i, false);
            assert!(self.is_consistent_copy_num(&new_a));
            assert!(self.is_consistent_copy_num(&new_b));
        }
    }
    pub fn cycle_and_direction_candidates(&self, copy_nums: &[u32]) -> Vec<(usize, bool)> {
        let candidates: Vec<(usize, bool)> = iproduct!((0..self.n_cycles()), &[true, false])
            .filter(|(cycle_id, &is_up)| self.is_acceptable(copy_nums, *cycle_id, is_up))
            .map(|(cycle_id, is_up)| (cycle_id, *is_up))
            .collect();
        candidates
    }
    pub fn from_seqs(seqs: &[Vec<u8>], k: usize) -> (CompressedDBG, Vec<u32>) {
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
    pub fn true_copy_nums_from_seqs(&self, seqs: &[Vec<u8>], k: usize) -> Option<Vec<u32>> {
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
    /// Get sequences that fully represents this cdbg
    /// by walking Eulerian circuit
    pub fn to_seqs(&self, copy_nums: &[u32]) -> Vec<Vec<u8>> {
        // counts[node] = (remaining copy_num of the node(k-mer))
        let mut counts = copy_nums.to_vec();
        let mut seqs = Vec::new();

        // find nodes that is not visited yet
        for head in self.heads().iter() {
            if counts[head.0] > 0 {
                seqs.push(self.travarse_kmers(head, &mut counts));
            }
        }
        // find nodes that is not visited yet
        loop {
            match counts.iter().position(|&c| c > 0) {
                Some(i) => {
                    let node = Node(i);
                    seqs.push(self.travarse_kmers(&node, &mut counts));
                }
                None => break,
            }
        }

        // check count is all zero
        let is_zero = counts.iter().all(|&c| c == 0);
        assert!(is_zero);
        debug!("is_zero: {}", is_zero);
        seqs
    }
    /// Generates string representation of this cdbg
    pub fn to_seqs_string(&self, copy_nums: &[u32]) -> String {
        let seqs_string: Vec<String> = self
            .to_seqs(copy_nums)
            .iter()
            .map(|seq| seq.iter().map(|&s| s as char).collect::<String>())
            .collect();
        seqs_string.join(",")
    }
    fn travarse_kmers(&self, from: &Node, counts: &mut [u32]) -> Vec<u8> {
        debug!("start: {}", self.kmer(from));
        // add from (first kmer)
        let mut now = from;
        let mut seq = self.kmer(from).to_vec();
        counts[from.0] -= 1;
        loop {
            match self
                .childs(now)
                .iter()
                .filter(|child| counts[child.0] > 0)
                .next()
            {
                Some(child) => {
                    debug!("next: {}", self.kmer(child));
                    // append a child (preceding kmer)
                    seq.push(self.emission(child));
                    counts[child.0] -= 1;
                    now = child;
                }
                None => break,
            };
        }
        seq
    }
    pub fn check_kmer_existence(&self, seqs: &[Vec<u8>], k: usize) -> (u32, u32) {
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
        (t, f)
    }
    /// prior score of this
    /// Assuming genome size ~ Normal(ave, std)
    pub fn prior_score(&self, copy_nums: &[u32], ave_size: u32, std_size: u32) -> Prob {
        normal_bin(self.total_emitable_copy_num(copy_nums), ave_size, std_size)
    }

    /// freqs -> trans_probs conversion function
    /// useful to create PHMM
    pub fn freq_to_trans_prob(&self, freqs: &[f64]) -> Vec<Vec<Prob>> {
        assert_eq!(self.n_kmers(), freqs.len());
        self.iter_nodes()
            .map(|v| {
                let fs: Vec<f64> = self
                    .childs(&v)
                    .iter()
                    .map(|&w| {
                        if self.is_emitable(&w) {
                            freqs[w.0]
                        } else {
                            0.0
                        }
                    })
                    .collect();
                let total_f: f64 = fs.iter().sum();
                fs.iter()
                    .map(|&f| {
                        if self.is_emitable(&v) && total_f > 0.0 {
                            Prob::from_prob(f / total_f)
                        } else {
                            Prob::from_prob(0.0)
                        }
                    })
                    .collect()
            })
            .collect()
    }
    /// calc the sum of freqs
    pub fn total_emitable_freq(&self, freqs: &[f64]) -> f64 {
        (0..self.n_kmers())
            .map(|i| Node(i))
            .filter(|&v| self.is_emitable(&v))
            .map(|v| freqs[v.0])
            .sum()
    }
    pub fn copy_nums_to_freqs(&self, freqs: &[u32]) -> Vec<f64> {
        freqs.iter().map(|&f| f as f64).collect()
    }

    pub fn to_edge_centric_graph(&self) -> (HashMap<Kmer, Node>, Vec<(Node, Node)>) {
        // assign index to all prefix/suffix
        let mut nodes: HashMap<Kmer, Node> = HashMap::default();
        // register NNNNN as Node(0), because it is useful as a start point
        nodes.insert(null_kmer(self.k - 1), Node(0));
        // register other k-1 mers
        for v in self.iter_nodes() {
            let kmer = self.kmer(&v);
            // add prefix
            let v = Node(nodes.len());
            nodes.entry(kmer.prefix()).or_insert(v);

            // add suffix
            let v = Node(nodes.len());
            nodes.entry(kmer.suffix()).or_insert(v);
        }

        // create edge for each kmer i.e. kmer -> (prefix, suffix)
        let edges: Vec<(Node, Node)> = self
            .iter_nodes()
            .map(|v| {
                let kmer = self.kmer(&v);
                let prefix = *nodes.get(&kmer.prefix()).unwrap();
                let suffix = *nodes.get(&kmer.suffix()).unwrap();
                (prefix, suffix)
            })
            .collect();

        (nodes, edges)
    }
    /// collect prefix/suffix as node, and assign the index
    /// graph with forward and backward edge
    pub fn to_indexed_digraph(&self) -> IndexedDiGraph {
        let (_, edges_forward) = self.to_edge_centric_graph();

        // reverse edges will be added
        let edges_reverse: Vec<(Node, Node)> = edges_forward
            .iter()
            .map(|&(prefix, suffix)| (suffix, prefix))
            .collect();

        let edges: Vec<_> = edges_forward
            .into_iter()
            .zip(edges_reverse.into_iter())
            .map(|(ef, er)| vec![ef, er])
            .flatten()
            .collect();

        IndexedDiGraph::from(edges)
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
    /// dot with probability (score) on each node
    pub fn as_dot_with_probs(&self, probs: &[Prob]) -> String {
        let mut s = String::new();
        writeln!(&mut s, "digraph cdbg {{");
        for v in self.iter_nodes() {
            // for node
            let kmer = self.kmer(&v);
            let prob = probs[v.0];
            writeln!(&mut s, "\t{} [label=\"{} {}\"];", v.0, kmer, prob);
            // for edges
            for w in self.childs(&v).iter() {
                writeln!(&mut s, "\t{} -> {};", v.0, w.0);
            }
        }
        writeln!(&mut s, "}}");
        s
    }
    /// show a histogram of cycle length distribution
    pub fn as_cycle_histogram(&self) -> String {
        let mut s = String::new();
        writeln!(&mut s, "#cycles={}", self.n_cycles());
        let mut histogram = Histogram::with_buckets(10);
        for i in 0..self.n_cycles() {
            histogram.add(self.cycle_components(i).len() as u64);
        }
        writeln!(&mut s, "{}", histogram);
        s
    }
    // stats related
    pub fn as_dbg_stats(&self) -> stats::DbgStats {
        let mut r = stats::DbgStats::default();
        r.k = self.k;
        r.n_kmers = self.n_kmers as u32;
        r.n_starting_kmers = self.starting_kmers().len() as u32;
        r
    }
    pub fn as_degree_stats(&self) -> stats::DegreeStats {
        let mut r = stats::DegreeStats::default();
        for v in self.iter_nodes() {
            let in_deg = self.parents(&v).len();
            let out_deg = self.childs(&v).len();
            r.in_degs[in_deg] += 1;
            r.out_degs[out_deg] += 1;
        }
        r
    }
    pub fn as_copy_num_stats(&self, copy_nums: &[u32]) -> stats::CopyNumStats {
        let mut r = stats::CopyNumStats::default();
        r.is_consistent = self.is_consistent_copy_num(copy_nums);
        r.total_emitable = self.total_emitable_copy_num(copy_nums);
        r.min = *copy_nums.iter().min().unwrap();
        r.max = *copy_nums.iter().max().unwrap();
        r.total = copy_nums.iter().sum();
        r.average = r.total as f32 / copy_nums.len() as f32;
        r.n_zero_copy_kmer = copy_nums.iter().filter(|c| **c == 0).count() as u32;
        r.n_nonzero_copy_kmer = copy_nums.len() as u32 - r.n_zero_copy_kmer;
        r
    }
    pub fn as_cycle_summary_stats(&self) -> stats::CycleSummaryStats {
        let mut r = stats::CycleSummaryStats::default();
        r.n_cycles = self.n_cycles() as u32;
        let cycle_lengths: Vec<u32> = (0..self.n_cycles())
            .map(|i| self.cycle_components(i).len() as u32)
            .collect();
        if self.n_cycles() == 0 {
            r.min_cycle_len = 0;
            r.max_cycle_len = 0;
            r.average_cycle_len = 0f32;
        } else {
            r.min_cycle_len = *cycle_lengths.iter().min().unwrap();
            r.max_cycle_len = *cycle_lengths.iter().max().unwrap();
            let total: u32 = cycle_lengths.iter().sum();
            r.average_cycle_len = total as f32 / self.n_cycles() as f32;
        }
        r
    }
    pub fn as_cycle_stats(&self, cycle_id: usize) -> stats::CycleStats {
        let mut r = stats::CycleStats::default();
        r.id = cycle_id;
        r.len = self.cycle_components(cycle_id).len();
        r.n_reverse = self
            .cycle_components(cycle_id)
            .iter()
            .map(|(_, dir)| match dir {
                cycles::CycleDirection::Reverse => 1,
                _ => 0,
            })
            .sum();
        r
    }
    pub fn as_all_stats(&self, copy_nums: &[u32]) -> stats::AllStats {
        let mut r = stats::AllStats::default();
        r.dbg = Some(self.as_dbg_stats());
        r.copy_num = Some(self.as_copy_num_stats(copy_nums));
        r.degree = Some(self.as_degree_stats());
        r.cycle_summary = Some(self.as_cycle_summary_stats());
        let cycles = (0..self.n_cycles())
            .map(|i| self.as_cycle_stats(i))
            .collect();
        r.cycles = Some(cycles);
        r
    }
    pub fn as_true_kmer_stats(
        &self,
        copy_nums: &[u32],
        copy_nums_true: &[u32],
    ) -> stats::TrueKmerStats {
        let mut r = stats::TrueKmerStats::default();
        r.size = self.total_emitable_copy_num(&copy_nums);
        r.true_size = self.total_emitable_copy_num(&copy_nums_true);
        let true_kmers: Vec<(u32, u32)> = copy_nums
            .iter()
            .zip(copy_nums_true.iter())
            .filter(|(&cn, &cnt)| cn > 0 && cnt > 0)
            .map(|(&cn, &cnt)| (cn, cnt))
            .collect();
        let false_kmers: Vec<(u32, u32)> = copy_nums
            .iter()
            .zip(copy_nums_true.iter())
            .filter(|(&cn, &cnt)| cn > 0 && cnt == 0)
            .map(|(&cn, &cnt)| (cn, cnt))
            .collect();
        r.n_true_kmer = true_kmers.iter().count();
        r.n_false_kmer = false_kmers.iter().count();

        r.true_kmer_max_copy_num = *true_kmers.iter().map(|(cn, _)| cn).max().unwrap_or(&0);
        r.true_kmer_min_copy_num = *true_kmers.iter().map(|(cn, _)| cn).min().unwrap_or(&0);
        r.true_kmer_ave_copy_num =
            true_kmers.iter().map(|(cn, _)| cn).sum::<u32>() as f32 / r.n_true_kmer as f32;

        r.false_kmer_max_copy_num = *false_kmers.iter().map(|(cn, _)| cn).max().unwrap_or(&0);
        r.false_kmer_min_copy_num = *false_kmers.iter().map(|(cn, _)| cn).min().unwrap_or(&0);
        r.false_kmer_ave_copy_num =
            false_kmers.iter().map(|(cn, _)| cn).sum::<u32>() as f32 / r.n_false_kmer as f32;

        r.copy_nums = copy_nums.to_vec();
        r.copy_nums_true = copy_nums_true.to_vec();

        r
    }
    pub fn layout_2d(&self) -> Vec<(Kmer, Pos)> {
        let idg = self.to_indexed_digraph();
        let positions = idg.layout_by_force_atlas2();
        let is_used: Vec<bool> = vec![false; idg.n_nodes()];
        let mut layout = Vec::new();
        for v in self.iter_nodes() {
            let kmer = self.kmer(&v);
            let (v, w) = idg.node_pair(&Edge(v.0));
            // add v
            if !is_used[v.0] {
                layout.push((kmer.prefix(), positions[v.0]))
            }
            // add w
            if !is_used[w.0] {
                layout.push((kmer.suffix(), positions[w.0]))
            }
        }
        layout
    }
    pub fn dump_layout_2d(&self) {
        for (km1mer, pos) in self.layout_2d().iter() {
            println!("N\t{}\t{}\t{}", km1mer, pos.0, pos.1);
        }
    }
    pub fn dump_edge_list(&self, copy_nums: &[u32]) {
        for v in self.iter_nodes() {
            println!("E\t{}\t{}", self.kmer(&v), copy_nums[v.0]);
        }
    }
    /// Create json to feed cytoscape.js graph
    /// format is array of either:
    ///   { group: 'nodes', data: { id } }
    /// or
    ///   { group: 'edges', data: { id, source, target, widths } }
    pub fn to_cytoscape_elements(
        &self,
        copy_nums_list: &[Vec<u32>],
        true_copy_nums: Option<&[u32]>,
    ) -> Vec<Element> {
        let (nodes, edges) = self.to_edge_centric_graph();
        let mut elements = Vec::new();

        // nodes
        let n_nodes = nodes.len();
        for (kmer, node) in nodes.into_iter() {
            elements.push(Element::Node {
                id: node.0,
                label: kmer,
            });
        }

        // edges
        for (i, &(v, w)) in edges.iter().enumerate() {
            let true_width = match true_copy_nums {
                Some(c) => Some(c[i]),
                None => None,
            };
            elements.push(Element::Edge {
                id: i + n_nodes, // element id should be unique, over both nodes and edges
                source: v.0,
                target: w.0,
                label: self.kmer(&Node(i)).clone(),
                widths: copy_nums_list
                    .iter()
                    .map(|copy_nums| copy_nums[i])
                    .collect(),
                true_width,
            });
        }

        elements
    }
    pub fn to_cytoscape_json_with_true(
        &self,
        copy_nums_list: &[Vec<u32>],
        true_copy_nums: Option<&[u32]>,
    ) -> String {
        let elements = self.to_cytoscape_elements(copy_nums_list, true_copy_nums);
        // serde_json::to_string_pretty(&elements).unwrap()
        serde_json::to_string(&elements).unwrap()
    }
    pub fn to_cytoscape_json(&self, copy_nums_list: &[Vec<u32>]) -> String {
        self.to_cytoscape_json_with_true(copy_nums_list, None)
    }
    pub fn to_gexf(&self) -> String {
        let mut s = String::new();
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn starting() {
        let mut d = DbgHash::new();
        d.add_seq(b"ATCGATTCGATCGATTCGATAGATCG", 8);
        let cdbg = CompressedDBG::from(&d, 8);
        let starting = cdbg.starting_kmers();
        // only 1 starting kmer
        assert_eq!(starting.len(), 1);
        // starting kmer is the first 7-mer
        assert_eq!(cdbg.kmer(&starting[0]).clone(), Kmer::from(b"NATCGATT"));
    }

    #[test]
    fn equality() {
        // same
        let mut d = DbgHash::new();
        d.add_seq(b"ATCGATTCGATCGATTCGATAGATCG", 8);
        let cdbg1a = CompressedDBG::from(&d, 8);
        let cdbg1b = CompressedDBG::from(&d, 8);
        assert_eq!(cdbg1a.as_dot(), cdbg1b.as_dot());
        assert_eq!(cdbg1a, cdbg1b);

        // different
        let mut d2 = DbgHash::new();
        d2.add_seq(b"ATCGATTCTTTAGTATTCGATAGATCG", 8);
        let cdbg2 = CompressedDBG::from(&d2, 8);
        assert_ne!(cdbg1a.as_dot(), cdbg2.as_dot());
        assert_ne!(cdbg1a, cdbg2);
    }

    #[test]
    fn cycle_candidates() {
        // sample data, loop less dbg
        let seqs = vec![b"ATTCGATCGATTT".to_vec()];
        let (cdbg, cn) = CompressedDBG::from_seqs(&seqs, 8);
        let cn_double: Vec<u32> = cn.iter().map(|x| x * 2).collect();
        let cn_zero: Vec<u32> = cn.iter().map(|_| 0).collect();

        assert_eq!(cdbg.n_cycles(), 1);
        assert_eq!(cdbg.cycle_and_direction_candidates(&cn).len(), 2);
        assert_eq!(cdbg.cycle_and_direction_candidates(&cn_double).len(), 2);
        assert_eq!(cdbg.cycle_and_direction_candidates(&cn_zero).len(), 1);
    }

    #[test]
    fn to_indexed_digraph() {
        let seqs = vec![b"ATTCGATCGATTT".to_vec()];
        let (cdbg, cn) = CompressedDBG::from_seqs(&seqs, 8);
        let idg = cdbg.to_indexed_digraph();
        assert_eq!(cdbg.n_kmers() * 2, idg.n_edges());
    }

    #[test]
    fn to_cytoscape_json() {
        let seqs = vec![b"ATTCGAC".to_vec()];
        let (cdbg, cn) = CompressedDBG::from_seqs(&seqs, 3);
        assert_eq!(
            cdbg.to_cytoscape_json(&[]),
            r#"[{"group":"nodes","data":{"id":7,"label":"NA"}},{"group":"nodes","data":{"id":6,"label":"CG"}},{"group":"nodes","data":{"id":5,"label":"AT"}},{"group":"nodes","data":{"id":3,"label":"AC"}},{"group":"nodes","data":{"id":8,"label":"GA"}},{"group":"nodes","data":{"id":4,"label":"CN"}},{"group":"nodes","data":{"id":1,"label":"TT"}},{"group":"nodes","data":{"id":2,"label":"TC"}},{"group":"nodes","data":{"id":0,"label":"NN"}},{"group":"edges","data":{"id":9,"source":1,"target":2,"label":"TTC","widths":[],"true_width":null}},{"group":"edges","data":{"id":10,"source":3,"target":4,"label":"ACN","widths":[],"true_width":null}},{"group":"edges","data":{"id":11,"source":5,"target":1,"label":"ATT","widths":[],"true_width":null}},{"group":"edges","data":{"id":12,"source":4,"target":0,"label":"CNN","widths":[],"true_width":null}},{"group":"edges","data":{"id":13,"source":2,"target":6,"label":"TCG","widths":[],"true_width":null}},{"group":"edges","data":{"id":14,"source":0,"target":7,"label":"NNA","widths":[],"true_width":null}},{"group":"edges","data":{"id":15,"source":6,"target":8,"label":"CGA","widths":[],"true_width":null}},{"group":"edges","data":{"id":16,"source":7,"target":5,"label":"NAT","widths":[],"true_width":null}},{"group":"edges","data":{"id":17,"source":8,"target":3,"label":"GAC","widths":[],"true_width":null}}]"#
        )
    }
}
