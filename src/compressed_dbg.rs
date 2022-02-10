use crate::cycles;
use crate::cycles::CycleDirection;
use crate::dbg::{DbgHash, DBG};
use crate::distribution::normal_bin;
use crate::graph::{IndexedDiGraph, Node};
use crate::kmer::kmer::{linear_seq_to_kmers, null_kmer, Kmer, KmerLike};
use crate::prob::Prob;
use fnv::FnvHashMap as HashMap;
use log::{debug, warn};
pub mod cycle;
pub mod output;

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
                    Some(&_v) => t += 1,
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
