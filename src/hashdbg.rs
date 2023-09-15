//!
//! HashDbg
//!
use crate::common::{CopyNum, Reads, Seq, StyledSequence};
use crate::graph::{
    compact::compact_simple_paths_for_targeted_nodes,
    utils::{degree_stats, split_node},
};
use crate::kmer::kmer::{
    linear_fragment_sequence_to_kmers, sequence_to_kmers, styled_sequence_to_kmers, Kmer, KmerLike,
};
use crate::multi_dbg::draft::{ErrorMetric, MinSquaredErrorCopyNumAndFreq, V1Error};
use fnv::{FnvHashMap as HashMap, FnvHashSet as HashSet};
use itertools::Itertools;
use petgraph::{
    algo::kosaraju_scc,
    graph::{DiGraph, EdgeIndex, NodeIndex},
    visit::EdgeRef,
    Direction,
};
use rustflow::min_flow::min_cost_flow_convex_fast;
use std::iter::Iterator;

///
/// De Bruijn graph structure that is implemented with a HashMap storing
/// `KmerLike -> CopyNum` mapping.
///
#[derive(Debug)]
pub struct HashDbg<K: KmerLike> {
    k: usize,
    kmers: HashMap<K, CopyNum>,
}

///
/// Basic operations
///
impl<K: KmerLike> HashDbg<K> {
    /// Create new HashDbg with no k-mers
    pub fn new(k: usize) -> HashDbg<K> {
        HashDbg {
            k,
            kmers: HashMap::default(),
        }
    }
    /// Create HashDbg from a set of k-mers with copy number
    pub fn from_kmers(k: usize, kmers_and_copy_nums: Vec<(K, CopyNum)>) -> HashDbg<K> {
        // assert all kmers are of the same size
        assert!(kmers_and_copy_nums.iter().all(|(kmer, _)| kmer.k() == k));

        // convert Vec into HashMap
        let kmers = kmers_and_copy_nums.into_iter().collect();
        HashDbg { k, kmers }
    }
    /// Size of k-mer in HashDbg
    pub fn k(&self) -> usize {
        self.k
    }
    /// Set copy number of k-mer (edge)
    pub fn set(&mut self, kmer: K, copy_num: CopyNum) {
        assert_eq!(kmer.k(), self.k());
        self.kmers.insert(kmer, copy_num);
    }
    /// Get copy number of k-mer (edge)
    /// if the k-mer does not exist in DBG, the copy number is zero.
    pub fn get(&self, kmer: &K) -> CopyNum {
        assert_eq!(kmer.k(), self.k());
        match self.kmers.get(&kmer) {
            Some(&copy_num) => copy_num,
            _ => 0,
        }
    }
    /// Add k-mer count
    pub fn add(&mut self, kmer: K, copy_num: CopyNum) {
        let copy_num_old = self.get(&kmer);
        self.set(kmer, copy_num + copy_num_old);
    }
    /// Remove k-mer from DBG
    pub fn remove(&mut self, kmer: &K) {
        self.kmers.remove(kmer);
    }
    /// Check k-mer (edge) exists in DBG
    pub fn has(&self, kmer: &K) -> bool {
        assert_eq!(kmer.k(), self.k());
        self.kmers.contains_key(kmer)
    }
    /// Get a list of edges (k-mers)
    pub fn edges(&self) -> Vec<K> {
        self.kmers.keys().cloned().collect()
    }
    /// Get a list of nodes (k-1-mers)
    pub fn nodes(&self) -> Vec<K> {
        let mut nodes = HashSet::default();
        for kmer in self.kmers.keys() {
            nodes.insert(kmer.prefix());
            nodes.insert(kmer.suffix());
        }
        nodes.into_iter().collect()
    }
    /// Get the number of k-mers (edges) in the DBG
    pub fn n(&self) -> usize {
        self.kmers.len()
    }
    pub fn kmers(&self) -> impl Iterator<Item = &K> {
        self.kmers.keys()
    }
    /// Source node (k-1-mer) from edge (k-mer)
    pub fn source(&self, kmer: &K) -> K {
        assert!(self.has(kmer));
        kmer.prefix()
    }
    /// Target node (k-1-mer) from edge (k-mer)
    pub fn target(&self, kmer: &K) -> K {
        assert!(self.has(kmer));
        kmer.suffix()
    }
    /// List of incoming edges of node (k-1-mer)
    /// km1mer XX -> kmers [YXX] whose suffix is the km1mer
    pub fn edges_in(&self, km1mer: &K) -> Vec<K> {
        assert_eq!(km1mer.k(), self.k() - 1);
        km1mer
            .preds()
            .into_iter()
            .filter(|kmer| self.has(kmer))
            .collect()
    }
    /// List of outgoing edges of node (k-1-mer)
    /// km1mer XX -> kmers [XXY] whose prefix is the km1mer
    pub fn edges_out(&self, km1mer: &K) -> Vec<K> {
        assert_eq!(km1mer.k(), self.k() - 1);
        km1mer
            .succs()
            .into_iter()
            .filter(|kmer| self.has(kmer))
            .collect()
    }
    /// Deprecated
    pub fn childs(&self, kmer: &K) -> Vec<K> {
        assert_eq!(kmer.k(), self.k());
        kmer.childs()
            .into_iter()
            .filter(|child| self.has(child))
            .collect()
    }
    /// Deprecated
    pub fn parents(&self, kmer: &K) -> Vec<K> {
        assert_eq!(kmer.k(), self.k());
        kmer.parents()
            .into_iter()
            .filter(|parent| self.has(parent))
            .collect()
    }
    /// Deprecated
    pub fn siblings(&self, kmer: &K) -> Vec<K> {
        assert_eq!(kmer.k(), self.k());
        kmer.siblings()
            .into_iter()
            .filter(|sibling| self.has(sibling))
            .collect()
    }
    ///
    /// add all kmers in linear seq (with leading/trailing NNN kmers)
    ///
    pub fn add_seq(&mut self, seq: &[u8]) {
        for kmer in sequence_to_kmers(seq, self.k()) {
            self.add(kmer, 1);
        }
    }
    ///
    /// add all kmers in linear fragment seq (without NNN kmers)
    ///
    pub fn add_fragment_seq(&mut self, seq: &[u8]) {
        for kmer in linear_fragment_sequence_to_kmers(seq, self.k()) {
            self.add(kmer, 1);
        }
    }
    ///
    /// add all kmers in styled sequence
    ///
    pub fn add_styled_sequence(&mut self, s: &StyledSequence) {
        for kmer in styled_sequence_to_kmers(s, self.k()) {
            self.add(kmer, 1);
        }
    }
}

//
// Display
//
impl<K: KmerLike + std::fmt::Display> std::fmt::Display for HashDbg<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // iter returns reference
        for kmer in self.kmers() {
            writeln!(f, "{} {}", kmer, self.get(&kmer))?;
        }
        Ok(())
    }
}

///
/// Constructors
///
/// * new
/// * from_kmers
///
/// * from_styled_seqs (like genome)
/// * from_fragment_seqs (like reads)
///
impl<K: KmerLike> HashDbg<K> {
    pub fn from_seq<S: Seq>(k: usize, seq: &S) -> Self {
        let mut d = HashDbg::new(k);
        d.add_seq(seq.as_ref());
        d
    }
    pub fn from_seqs<T>(k: usize, seqs: T) -> Self
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        let mut d = HashDbg::new(k);
        for seq in seqs {
            let seq = seq.as_ref();
            // ignore read if it is shorter than k
            if seq.len() >= k {
                d.add_seq(seq);
            }
        }
        d
    }
    /// Construct HashDbg from seqs regarding as Fragments
    /// Ignoring heads and tails (like nnnA and Gnnn).
    ///
    pub fn from_fragment_seqs<T>(k: usize, seqs: T) -> Self
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        let mut d = HashDbg::new(k);
        for seq in seqs {
            let seq = seq.as_ref();
            if seq.len() >= k {
                d.add_fragment_seq(seq);
            }
        }
        d
    }
    pub fn from_styled_seqs<T>(k: usize, seqs: T) -> Self
    where
        T: IntoIterator,
        T::Item: AsRef<StyledSequence>,
    {
        let mut d = HashDbg::new(k);
        for seq in seqs {
            d.add_styled_sequence(seq.as_ref());
        }
        d
    }
}

impl<K: KmerLike> HashDbg<K> {
    ///
    /// to kmer count profile
    ///
    pub fn to_kmer_profile(&self) -> HashMap<K, CopyNum> {
        self.kmers.clone()
    }
    /// Generate degree statistics of the graph
    /// as `HashMap<(in_degree: usize, out_degree), count_of_node: usize>`
    pub fn degree_stats(&self) -> HashMap<(usize, usize), usize> {
        let graph = self.to_graph(|_| (), |_| ());
        degree_stats(&graph)
    }
    /// Generate k-mer copy number summary
    /// as `HashMap<copy_num: usize, count: usize>`
    pub fn copy_num_stats(&self) -> HashMap<usize, usize> {
        let mut h = HashMap::default();
        for (_kmer, &copy_num) in self.kmers.iter() {
            *h.entry(copy_num).or_insert(0) += 1;
        }
        h
    }
    ///
    /// to edge-centric (full) DBG as petgraph::DiGraph
    /// edge is k-mer and node is k-1-mer.
    ///
    /// to_node(km1mer: K) -> N
    /// to_edge(kmer: K) -> E
    /// gives petgraph::DiGraph<N, E>
    ///
    pub fn to_graph<N, E, FN, FE>(&self, to_node: FN, to_edge: FE) -> DiGraph<N, E>
    where
        FN: Fn(&K) -> N,
        FE: Fn(&K) -> E,
    {
        let mut graph = DiGraph::new();

        // mapping from km1mer to node index
        let mut ids: HashMap<K, NodeIndex> = HashMap::default();

        // add each node (km1mer)
        for node in self.nodes() {
            let id = graph.add_node(to_node(&node));
            ids.insert(node, id);
        }

        // add each edge (kmer) between source/target nodes (km1mer)
        for edge in self.edges() {
            let s = *ids.get(&self.source(&edge)).unwrap();
            let t = *ids.get(&self.target(&edge)).unwrap();
            graph.add_edge(s, t, to_edge(&edge));
        }

        graph
    }
    /// Check the copy numbers are consistent
    /// i.e. for all nodes, sum of copy numbers of in-edges/out-edges is the same.
    pub fn is_copy_nums_consistent(&self) -> bool {
        self.nodes().into_iter().all(|node| {
            let copy_nums_in: CopyNum = self
                .edges_in(&node)
                .into_iter()
                .map(|edge_in| self.get(&edge_in))
                .sum();
            let copy_nums_out: CopyNum = self
                .edges_out(&node)
                .into_iter()
                .map(|edge_out| self.get(&edge_out))
                .sum();
            copy_nums_in == copy_nums_out
        })
    }
    /// Remove kmers whose count is less than `min_copy_num`
    /// Returns the number of removed k-mers
    pub fn remove_rare_kmers(&mut self, min_copy_num: CopyNum) -> usize {
        let n0 = self.n();
        self.kmers
            .retain(|_kmer, copy_num| *copy_num >= min_copy_num);
        let n1 = self.n();
        n0 - n1
    }
    /// k-mer (edge) is deadend (= edge that has no child/parent edges) or not?
    pub fn is_deadend(&self, kmer: &K) -> bool {
        self.childs(kmer).len() == 0 || self.parents(kmer).len() == 0
    }
    /// Remove deadend k-mers whose count is less than `min_count`
    ///
    /// Removing a deadend can create some new deadends that is childs/parents of that node
    ///
    /// Returns the number of removed deadends
    pub fn remove_deadends(&mut self, min_count: CopyNum) -> usize {
        let mut deadends: Vec<K> = self
            .edges()
            .into_iter()
            .filter(|edge| self.get(edge) < min_count)
            .filter(|edge| self.is_deadend(edge))
            .collect();
        let mut n_deadends = 0;
        println!("initial deadends {}", deadends.len());

        while let Some(deadend) = deadends.pop() {
            // remove the deadend
            self.remove(&deadend);
            n_deadends += 1;

            // child/parent of the deadend is added to the list if it is a new deadend
            for child in self.childs(&deadend) {
                if self.is_deadend(&child) {
                    deadends.push(child);
                }
            }
            for parent in self.parents(&deadend) {
                if self.is_deadend(&parent) {
                    deadends.push(parent);
                }
            }
        }

        println!("removed {} deadends", n_deadends);
        n_deadends
    }
    /// Add all starting k-mers of the k-mer (i.e. connect from terminal `nnn` into the k-mer).
    ///
    /// Example: if kmer is `AGCT`, then `nnnA, nnAG, nAGC` will be added.
    pub fn add_starting_kmers(&mut self, kmer: &K) {
        let copy_num = self.get(kmer);
        for starting_kmer in kmer.starting_kmers() {
            self.add(starting_kmer, copy_num);
        }
    }
    /// Add all ending k-mers of the k-mer (i.e. connect from the k-mer into terminal `nnn`).
    ///
    /// Example: if kmer is `AGCT`, then `GCTn, CTnn, Tnnn` will be added.
    pub fn add_ending_kmers(&mut self, kmer: &K) {
        let copy_num = self.get(kmer);
        for ending_kmer in kmer.ending_kmers() {
            self.add(ending_kmer, copy_num);
        }
    }
    /// Add starting/ending k-mers for all deadend k-mers.
    ///
    /// if k-mer has no parent, then add starting k-mers (nnn -> k-mer)
    /// if k-mer has no child, then add ending k-mers (k-mer -> nnn)
    ///
    /// Returns the starting/ending k-mers
    pub fn augment_deadends(&mut self) -> (Vec<K>, Vec<K>) {
        let mut starts = Vec::new();
        let mut ends = Vec::new();

        for edge in self.edges() {
            // (1) edge is starting point if it has no parent.
            // then add starting k-mers to the edge.
            if self.parents(&edge).len() == 0 {
                self.add_starting_kmers(&edge);
                starts.push(edge.clone());
            }

            // (2) edge is ending point if it has no child.
            // then add ending k-mers to the edge.
            if self.childs(&edge).len() == 0 {
                self.add_ending_kmers(&edge);
                ends.push(edge);
            }
        }

        (starts, ends)
    }
    /// DBG into connected components as a set of sets of k-mers.
    /// vec of components to be returned is sorted in descending order of its size.
    ///
    /// Implementation: convert HashDbg into petgraph::DiGraph and run `petgraph::algo::kosaraju_scc`
    pub fn connected_components(&self) -> Vec<Vec<K>> {
        // convert to petgraph::DiGraph whose edge weight is kmer: K
        let graph = self.to_graph(|_km1mer| (), |kmer| kmer.clone());
        // tarjan_scc
        let mut components = kosaraju_scc(&graph);
        components.sort_by_key(|component| std::cmp::Reverse(component.len()));
        // convert nodes in a component into in-coming edges
        // connected component defined by a set of nodes corresponds to a set of edges.
        components
            .into_iter()
            .map(|component| {
                component
                    .into_iter()
                    .flat_map(|node| {
                        graph
                            .edges_directed(node, Direction::Incoming)
                            .map(|edge| edge.weight().clone())
                    })
                    .collect()
            })
            .collect()
    }
    /// Get subgraph of this HashDbg
    pub fn largest_component(&self) -> HashDbg<K> {
        let mut components = self.connected_components();
        if components.is_empty() {
            // empty
            HashDbg::new(self.k())
        } else {
            let component = components.remove(0);
            let kmers: Vec<(K, CopyNum)> = component
                .into_iter()
                .map(|kmer| {
                    let copy_num = self.get(&kmer);
                    (kmer, copy_num)
                })
                .collect();
            HashDbg::from_kmers(self.k(), kmers)
        }
    }
    ///
    /// Find approximate copy numbers from k-mer counts
    ///
    /// * coverage
    /// * specify n_haplotypes
    /// * convert to min-flow network
    ///
    /// Return
    /// * Network flow graph (edge weight is `MinSquaredErrorCopyNumAndFreq`)
    /// * Mapping from edge index into k-mers as `Vec<Vec<K>>`
    ///
    pub fn to_min_squared_error_copy_nums_network<T: ErrorMetric>(
        &self,
        coverage: f64,
        n_haplotypes: Option<usize>,
    ) -> (DiGraph<(), MinSquaredErrorCopyNumAndFreq<T>>, Vec<Vec<K>>) {
        // convert to full graph
        // node weight is km1mer: K
        // edge weight is (kmer: K, count: CopyNum)
        let full = self.to_graph(
            |km1mer| km1mer.clone(),
            |kmer| (kmer.clone(), self.get(kmer)),
        );

        // compact
        // contract nodes except terminal node `nnn`
        let compact = compact_simple_paths_for_targeted_nodes(&full, |km1mer| !km1mer.is_null());

        // flow network
        let mut network = compact.map(
            |_node, _| (),
            |_edge, edge_weight| {
                let freqs = edge_weight
                    .iter()
                    .filter_map(|(_, (kmer, count))| {
                        if kmer.has_null() {
                            None
                        } else {
                            Some(*count as f64 / coverage)
                        }
                    })
                    .collect();
                MinSquaredErrorCopyNumAndFreq::new(freqs, None, false)
            },
        );

        // split terminal node into two to specify (fix) the number of haplotypes
        let terminal = compact
            .node_indices()
            .find(|&node| compact[node].is_null())
            .expect("graph has no terminal node");
        split_node(
            &mut network,
            terminal,
            Some(MinSquaredErrorCopyNumAndFreq::new(
                vec![],
                n_haplotypes,
                false,
            )),
        );

        // mapping from edge into kmers
        let mut kmers: Vec<Vec<K>> = compact
            .edge_indices()
            .map(|edge| {
                compact[edge]
                    .iter()
                    .map(|(_, (kmer, _))| kmer.clone())
                    .collect()
            })
            .collect();
        // the last edge added by splitting terminal does not correspond to any k-mer
        kmers.push(vec![]);

        (network, kmers)
    }
    ///
    pub fn generate_hashdbg_with_min_squared_error_copy_nums(
        &self,
        coverage: f64,
        n_haplotypes: Option<usize>,
    ) -> HashDbg<K> {
        println!("converting");
        let (network, kmer_map) =
            self.to_min_squared_error_copy_nums_network::<V1Error>(coverage, n_haplotypes);
        println!(
            "network n_nodes={} n_edges={}",
            network.node_count(),
            network.edge_count()
        );
        let copy_nums =
            min_cost_flow_convex_fast(&network).expect("mse flownetwork cannot be solved");

        let kmers: Vec<(K, CopyNum)> = network
            .edge_indices()
            .flat_map(|edge| {
                let copy_num = copy_nums[edge];
                kmer_map[edge.index()]
                    .iter()
                    .map(move |kmer| (kmer.clone(), copy_num))
            })
            .collect();
        HashDbg::from_kmers(self.k(), kmers)
    }
}

///
/// Output
///
impl<K: KmerLike> HashDbg<K> {
    ///
    ///
    pub fn to_gfa_writer<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        let graph = self.to_graph(
            |km1mer| km1mer.clone(),
            |kmer| (kmer.clone(), self.get(kmer)),
        );
        let compact = compact_simple_paths_for_targeted_nodes(&graph, |km1mer| !km1mer.is_null());

        for edge in compact.edge_indices() {
            let weight = &compact[edge];
            let label = weight.iter().map(|(_, (_, count))| count).join(",");
            let count_sum: usize = weight.iter().map(|(_, (_, count))| count).sum();
            let count_ave = count_sum as f64 / weight.len() as f64;
            // let seq = &self.seq_compact(edge);
            writeln!(
                writer,
                "S\t{}\t*\tDP:f:{}\tLN:i:{}\tLB:Z:{}",
                edge.index(),
                // sequence_to_string(&seq),
                count_ave,
                weight.len(),
                label,
            )?
        }
        for node in compact.node_indices() {
            let km1mer = &compact[node];
            if !km1mer.is_null() {
                for in_edge in compact.edges_directed(node, Direction::Incoming) {
                    for out_edge in compact.edges_directed(node, Direction::Outgoing) {
                        writeln!(
                            writer,
                            "L\t{}\t+\t{}\t+\t*\tID:Z:{}",
                            in_edge.id().index(),
                            out_edge.id().index(),
                            node.index(),
                        )?
                    }
                }
            }
        }
        Ok(())
    }
    ///
    ///
    pub fn to_gfa_file<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        let mut file = std::fs::File::create(path).unwrap();
        self.to_gfa_writer(&mut file)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::veckmer::{kmer, VecKmer};

    #[test]
    fn hashdbg_v2_new() {
        let mut hd: HashDbg<VecKmer> = HashDbg::new(4);

        // add get
        assert_eq!(hd.get(&kmer(b"ATCG")), 0);
        hd.add(kmer(b"ATCG"), 1);
        hd.add(kmer(b"TTTT"), 2);
        assert_eq!(hd.get(&kmer(b"TTTT")), 2);
        hd.add(kmer(b"TTTT"), 1);
        assert_eq!(hd.get(&kmer(b"TTTT")), 3);
        assert_eq!(hd.get(&kmer(b"ATCG")), 1);
        println!("{:?}", hd);

        // is_exists
        assert_eq!(hd.has(&kmer(b"ATCG")), true);
        assert_eq!(hd.has(&kmer(b"ATCT")), false);

        // kmers
        let kmers: Vec<VecKmer> = hd.kmers().cloned().collect();
        assert_eq!(kmers, vec![kmer(b"TTTT"), kmer(b"ATCG")]);
        assert_eq!(hd.edges(), vec![kmer(b"TTTT"), kmer(b"ATCG")]);
        assert_eq!(hd.nodes(), vec![kmer(b"ATC"), kmer(b"TCG"), kmer(b"TTT")]);
        assert_eq!(hd.edges_out(&kmer(b"ATC")), vec![kmer(b"ATCG")]);
        assert_eq!(hd.edges_in(&kmer(b"ATC")), vec![]);
        assert_eq!(hd.edges_out(&kmer(b"TTT")), vec![kmer(b"TTTT")]);
        assert_eq!(hd.edges_in(&kmer(b"TTT")), vec![kmer(b"TTTT")]);

        // set and delete
        hd.set(kmer(b"TTTT"), 0);
        assert_eq!(hd.has(&kmer(b"TTTT")), true);
        hd.remove(&kmer(b"TTTT"));
        assert_eq!(hd.has(&kmer(b"TTTT")), false);

        assert_eq!(hd.edges(), vec![kmer(b"ATCG")]);
        assert_eq!(hd.nodes(), vec![kmer(b"ATC"), kmer(b"TCG")]);
    }

    #[test]
    fn hashdbg_v2_seq() {
        let mut hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"ATCGATTCGAT");
        println!("{}", hd);

        // childs parents siblings
        assert_eq!(hd.childs(&kmer(b"ATCG")), vec![kmer(b"TCGA")]);
        assert_eq!(hd.childs(&kmer(b"ATCC")), vec![]);
        assert_eq!(hd.parents(&kmer(b"CGAT")), vec![kmer(b"TCGA")]);
        assert_eq!(
            hd.childs(&kmer(b"CGAT")),
            vec![kmer(b"GATT"), kmer(b"GATn")]
        );
        assert_eq!(hd.childs(&kmer(b"Tnnn")), vec![kmer(b"nnnA")]);
        assert_eq!(
            hd.siblings(&kmer(b"ATCG")),
            vec![kmer(b"ATCG"), kmer(b"TTCG")]
        );
        assert_eq!(hd.siblings(&kmer(b"nnnA")), vec![kmer(b"nnnA")]);

        assert_eq!(hd.edges().len(), 12);
        assert_eq!(hd.nodes().len(), 11);

        for e in hd.edges_out(&kmer(b"nnn")) {
            println!("out {}", e);
        }
        assert_eq!(hd.edges_out(&kmer(b"nnn")), vec![kmer(b"nnnA")]);
        assert_eq!(hd.edges_in(&kmer(b"nnn")), vec![kmer(b"Tnnn")]);
        assert_eq!(hd.edges_out(&kmer(b"TCG")), vec![kmer(b"TCGA")]);
        assert_eq!(
            hd.edges_out(&kmer(b"GAT")),
            vec![kmer(b"GATT"), kmer(b"GATn")]
        );
        assert_eq!(hd.get(&kmer(b"TCGA")), 2);
        assert_eq!(hd.source(&kmer(b"TCGA")), kmer(b"TCG"));
        assert_eq!(hd.target(&kmer(b"TCGA")), kmer(b"CGA"));

        //
        // graph conversion
        //
        let g = hd.to_graph(|km1mer| km1mer.clone(), |kmer| kmer.clone());
        assert_eq!(g.edge_count(), 12);
        assert_eq!(g.node_count(), 11);
        for e in g.edge_indices() {
            let (s, t) = g.edge_endpoints(e).unwrap();
            assert_eq!(g[s], g[e].prefix());
            assert_eq!(g[t], g[e].suffix());
        }
        println!("{}", petgraph::dot::Dot::with_config(&g, &[]));
    }

    #[test]
    fn hashdbg_v2_profile() {
        let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"ATCGATTCGAT");
        let p = hd.to_kmer_profile();
        for (k, v) in p.iter() {
            println!("{} {}", k, v);
        }

        // answer
        let mut p2 = HashMap::default();
        p2.insert(kmer(b"ATTC"), 1);
        p2.insert(kmer(b"TTCG"), 1);
        p2.insert(kmer(b"nnAT"), 1);
        p2.insert(kmer(b"CGAT"), 2);
        p2.insert(kmer(b"TCGA"), 2);
        p2.insert(kmer(b"nATC"), 1);
        p2.insert(kmer(b"GATT"), 1);
        p2.insert(kmer(b"GATn"), 1);
        p2.insert(kmer(b"nnnA"), 1);
        p2.insert(kmer(b"ATnn"), 1);
        p2.insert(kmer(b"ATCG"), 1);
        p2.insert(kmer(b"Tnnn"), 1);

        assert_eq!(p, p2);
    }
}
