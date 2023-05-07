//!
//! Output (serialization and deserialization) functions of MultiDbg
//!
//! # DBG format [`MultiDbg::to_dbg_writer`]
//!
//! ```text
//! K
//! N
//! E
//!
//! # paths format
//! P
//! C
//! ```
//!
//! # GFA format [`MultiDbg::to_gfa_writer`]
//!
//! ## without posterior
//!
//! ```text
//! # segment
//! #            copy_num
//! S 0 TTCTG    DP:f:1  LN:i:5
//! S 1 CCTAGCT  DP:f:1  LN:i:7
//!
//! # link
//! L 0 + 1 + *   ID:Z:1
//! ```
//!
//! ## with posterior
//!
//! ```text
//! # segment
//! #            copy_num          posterior_distribution_of_copy_num,true_copy_num
//! S 0 TTCTG    DP:f:1.1  LN:i:5  LB:Z:P(1x=0.001,2x=0.9),1x
//! S 1 CCTAGCT  DP:f:0.9  LN:i:7
//!
//! # link
//! L 0 + 1 + *   ID:Z:1
//!
//! # path
//! P
//! ```
//!
use super::posterior::Posterior;
use super::{CopyNums, MultiCompactEdge, MultiCompactNode, MultiDbg, MultiFullEdge, MultiFullNode};
use crate::common::{sequence_to_string, CopyNum, ReadCollection, Seq, NULL_BASE};
use crate::hmmv2::hint::{Mapping, Mappings};
use crate::prob::Prob;
use flate2::bufread::GzDecoder;
use flate2::write::GzEncoder;

use flate2::Compression;
use itertools::Itertools;
use petgraph::graph::{DefaultIx, DiGraph, EdgeIndex, NodeIndex};

///
/// debug print methods
///
impl MultiDbg {
    ///
    /// Dot file with each node/edge shown in Display serialization
    ///
    pub fn to_dot(&self) -> String {
        format!(
            "Full:\n{}Compact:\n{}",
            petgraph::dot::Dot::with_config(&self.graph_full(), &[]),
            petgraph::dot::Dot::with_config(&self.graph_compact(), &[]),
        )
    }
    ///
    /// Debug output each node/edge with kmer
    ///
    pub fn show_graph_with_kmer(&self) {
        println!("k={}", self.k());
        println!("genome_size={}", self.genome_size());
        println!("is_copy_nums_valid={}", self.is_copy_nums_valid());
        println!("degree_stats={:?}", self.degree_stats());
        println!("to_string={}", self.to_string());

        println!("Full:");
        println!("n_nodes={}", self.n_nodes_full());
        println!("n_edges={}", self.n_edges_full());
        for (node, weight) in self.nodes_full() {
            println!("v{}\t{}\t{}", node.index(), weight, self.km1mer_full(node));
        }
        for (edge, s, t, weight) in self.edges_full() {
            println!(
                "e{}\tv{}\tv{}\t{}\t{}",
                edge.index(),
                s.index(),
                t.index(),
                weight,
                self.kmer_full(edge)
            );
        }

        println!("Compact:");
        println!("n_nodes={}", self.n_nodes_compact());
        println!("n_edges={}", self.n_edges_compact());
        for (node, weight) in self.nodes_compact() {
            println!(
                "v{}\t{}\t{}",
                node.index(),
                weight,
                self.km1mer_compact(node)
            );
        }
        for (edge, s, t, weight) in self.edges_compact() {
            println!(
                "e{}\tv{}\tv{}\t{}\t{}",
                edge.index(),
                s.index(),
                t.index(),
                weight,
                self.kmer_compact(edge)
            );
        }
    }
}

///
/// DBG format
///
impl MultiDbg {
    ///
    /// create DBG string with `to_dbg_writer`
    ///
    pub fn to_dbg_string(&self) -> String {
        let mut writer = Vec::with_capacity(128);
        self.to_dbg_writer(&mut writer).unwrap();
        String::from_utf8(writer).unwrap()
    }
    ///
    /// create DBG file with `to_dbg_writer`
    ///
    pub fn to_dbg_file<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        let mut file = std::fs::File::create(path.as_ref()).unwrap();

        if path.as_ref().extension().is_some_and(|ext| ext == "gz") {
            let mut writer = GzEncoder::new(file, Compression::default());
            self.to_dbg_writer(&mut writer)?;
            writer.try_finish()
        } else {
            self.to_dbg_writer(&mut file)
        }
    }
    /// DBG format
    ///
    /// ```text
    /// # comment
    /// # K: kmer-size
    /// K 4
    ///
    /// # N: node
    /// # id km1mer
    /// N 0  NNN
    /// N 1  ATC
    ///
    /// # E: edge
    /// # id source target seq       copy_num edges
    /// E 0  0      1      ATCGATGCT 10       8,7,2,3,1
    /// E 0  0      1      ATCGATGCT 5        8,7,2,3,1
    ///
    /// # P: path (sequence of edges)
    /// P 0,5,2,3,1,2,4,2
    ///
    /// # C: copy numbers
    /// # copy_nums_vector
    /// C 1,1,1,0,0,0,0,0,1
    /// ```
    pub fn to_dbg_writer<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        writeln!(writer, "# {}", env!("GIT_HASH"))?;
        writeln!(writer, "K\t{}", self.k())?;
        for (node, weight) in self.nodes_compact() {
            writeln!(writer, "N\t{}\t{}", node.index(), self.km1mer_compact(node))?
        }
        for (edge, s, t, weight) in self.edges_compact() {
            writeln!(
                writer,
                "E\t{}\t{}\t{}\t{}\t{}\t{}",
                edge.index(),
                s.index(),
                t.index(),
                self.kmer_compact(edge),
                self.copy_num_of_edge_in_compact(edge),
                weight.edges_in_full().iter().map(|e| e.index()).format(","),
            )?
        }
        Ok(())
    }
    ///
    ///
    ///
    pub fn from_dbg_reader<R: std::io::BufRead>(reader: R) -> Self {
        let mut k = None;
        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        // sum of seq in compact graph and the number of edges in full graph
        let mut n_bases = 0;

        for line in reader.lines() {
            let text = line.unwrap();
            let first_char = text.chars().nth(0).unwrap();
            match first_char {
                'K' => {
                    let mut iter = text.split_whitespace();
                    iter.next().unwrap(); // 'K'
                    k = Some(iter.next().unwrap().parse().unwrap());
                }
                'N' => {
                    let mut iter = text.split_whitespace();
                    iter.next().unwrap(); // 'N'

                    let node: NodeIndex<DefaultIx> =
                        NodeIndex::new(iter.next().unwrap().parse().unwrap());
                    let km1mer = iter.next().unwrap().as_bytes().to_vec();

                    assert_eq!(nodes.len(), node.index(), "node is not sorted");
                    nodes.push((node, km1mer));
                }
                'E' => {
                    let k = k.unwrap();
                    let mut iter = text.split_whitespace();
                    iter.next().unwrap(); // 'E'

                    let edge: EdgeIndex<DefaultIx> =
                        EdgeIndex::new(iter.next().unwrap().parse().unwrap());
                    let s: NodeIndex<DefaultIx> =
                        NodeIndex::new(iter.next().unwrap().parse().unwrap());
                    let t: NodeIndex<DefaultIx> =
                        NodeIndex::new(iter.next().unwrap().parse().unwrap());

                    let mut kmer: Vec<u8> = iter.next().unwrap().as_bytes().to_vec();
                    let seq = kmer.split_off(k - 1);
                    n_bases += seq.len();

                    let copy_num: CopyNum = iter.next().unwrap().parse().unwrap();
                    let edges_in_full: Vec<EdgeIndex<DefaultIx>> = iter
                        .next()
                        .unwrap()
                        .split(',')
                        .map(|s| EdgeIndex::new(s.parse().unwrap()))
                        .collect();
                    assert_eq!(
                        edges_in_full.len(),
                        seq.len(),
                        "length of seq and edges_in_full is different"
                    );

                    assert_eq!(edges.len(), edge.index(), "edge is not sorted");
                    edges.push((edge, s, t, seq, copy_num, edges_in_full));
                }
                '#' => {} // pass
                _ => {}   // ignore
            }
        }

        // full
        let mut full = DiGraph::new();
        for (_, km1mer) in nodes.iter() {
            let is_terminal = km1mer.iter().all(|&x| x == NULL_BASE);
            full.add_node(MultiFullNode::new(is_terminal));
        }
        let mut edges_full = vec![None; n_bases];
        for (_, s, t, seq, copy_num, edges_in_full) in edges.iter() {
            // Compact:
            // s ------------------> t
            //    e0,e1,e2...
            //    ATT...
            //
            // into
            //
            // Full:
            //
            // s ----> v ----> v' ----> ... ----> v'' ----> t
            //   e0      e1       e2
            //   A       T        T
            //
            let n = seq.len();
            let mut w_prev = None;
            for i in 0..n {
                let base = seq[i];
                let edge_in_full = edges_in_full[i];

                let v = if i == 0 {
                    *s
                } else {
                    // previous w in i-1
                    w_prev.unwrap()
                };

                let w = if i == n - 1 {
                    *t
                } else {
                    // new node
                    full.add_node(MultiFullNode::new(false))
                };

                edges_full[edge_in_full.index()] =
                    Some((v, w, MultiFullEdge::new(base, *copy_num)));

                // w will be source (v) in i+1
                w_prev = Some(w);
            }
        }
        for (i, e) in edges_full.into_iter().enumerate() {
            let (source, target, weight) = e.expect("index of edge in full is wrong");
            let edge = full.add_edge(source, target, weight);
            assert_eq!(i, edge.index());
        }

        // compact
        let mut compact = DiGraph::new();
        for (_, km1mer) in nodes.iter() {
            let is_terminal = km1mer.iter().all(|&x| x == NULL_BASE);
            compact.add_node(MultiCompactNode::new(is_terminal));
        }
        for (_, s, t, _, _, edges_in_full) in edges.into_iter() {
            compact.add_edge(s, t, MultiCompactEdge::new(edges_in_full));
        }

        MultiDbg {
            k: k.expect("no K section"),
            full,
            compact,
        }
    }
    ///
    /// parse DBG string with `from_dbg_reader`
    ///
    pub fn from_dbg_str(s: &str) -> Self {
        Self::from_dbg_reader(s.as_bytes())
    }
    ///
    /// parse DBG file with `from_dbg_reader`
    ///
    pub fn from_dbg_file<P: AsRef<std::path::Path>>(path: P) -> Self {
        let file = std::fs::File::open(path.as_ref()).unwrap();
        let reader = std::io::BufReader::new(file);

        if path.as_ref().extension().is_some_and(|ext| ext == "gz") {
            let decoder = GzDecoder::new(reader);
            Self::from_dbg_reader(std::io::BufReader::new(decoder))
        } else {
            Self::from_dbg_reader(reader)
        }
    }
}

///
/// PATHS format
///
impl MultiDbg {
    ///
    /// create PATHS file with [`MultiDbg::to_paths_writer`]
    ///
    pub fn to_paths_file<F, PS, P>(&self, filename: F, paths: PS) -> std::io::Result<()>
    where
        F: AsRef<std::path::Path>,
        PS: AsRef<[P]>,
        P: AsRef<[EdgeIndex]>,
    {
        let mut file = std::fs::File::create(filename).unwrap();
        self.to_paths_writer(&mut file, paths)
    }
    ///
    /// create PATHS string with [`MultiDbg::to_paths_writer`]
    ///
    pub fn to_paths_string<PS, P>(&self, paths: PS) -> String
    where
        PS: AsRef<[P]>,
        P: AsRef<[EdgeIndex]>,
    {
        let mut writer = Vec::with_capacity(128);
        self.to_paths_writer(&mut writer, paths).unwrap();
        String::from_utf8(writer).unwrap()
    }
    ///
    /// create PATHS file
    ///
    pub fn to_paths_writer<W, PS, P>(&self, mut writer: W, paths: PS) -> std::io::Result<()>
    where
        W: std::io::Write,
        PS: AsRef<[P]>,
        P: AsRef<[EdgeIndex]>,
    {
        for path in paths.as_ref().into_iter() {
            writeln!(
                writer,
                "P\t{}",
                path.as_ref().into_iter().map(|e| e.index()).format(",")
            )?;
        }
        Ok(())
    }
    ///
    ///
    ///
    pub fn from_paths_reader<R: std::io::BufRead>(reader: R) -> Vec<Vec<EdgeIndex>> {
        let mut paths = Vec::new();

        for line in reader.lines() {
            let text = line.unwrap();
            let first_char = text.chars().nth(0).unwrap();
            if first_char == 'P' {
                let mut iter = text.split_whitespace();

                // first segment
                iter.next().unwrap(); // 'P'

                // second segment
                let path: Vec<EdgeIndex> = iter
                    .next()
                    .unwrap()
                    .split(',')
                    .map(|s| EdgeIndex::new(s.parse().unwrap()))
                    .collect();

                // TODO check that this is valid path of current dbg

                paths.push(path);
            }
        }

        paths
    }
    ///
    /// parse PATHS string with [`MultiDbg::from_paths_reader`]
    ///
    pub fn from_paths_str(s: &str) -> Vec<Vec<EdgeIndex>> {
        Self::from_paths_reader(s.as_bytes())
    }
    ///
    /// parse PATHS file with [`MultiDbg::from_paths_reader`]
    ///
    pub fn from_paths_file<P: AsRef<std::path::Path>>(path: P) -> Vec<Vec<EdgeIndex>> {
        let file = std::fs::File::open(path).unwrap();
        let reader = std::io::BufReader::new(file);
        Self::from_paths_reader(reader)
    }
}

///
/// Mapping file (.MAP format file) dump/load functions
///
impl MultiDbg {
    ///
    /// create MAP file with [`MultiDbg::to_map_writer`]
    ///
    pub fn to_map_file<F: AsRef<std::path::Path>, S: Seq>(
        &self,
        filename: F,
        reads: &ReadCollection<S>,
        mappings: &Mappings,
    ) -> std::io::Result<()> {
        let mut file = std::fs::File::create(filename).unwrap();
        self.to_map_writer(&mut file, reads, mappings)
    }
    ///
    /// create MAP string with [`MultiDbg::to_map_writer`]
    ///
    pub fn to_map_string<S: Seq>(&self, reads: &ReadCollection<S>, mappings: &Mappings) -> String {
        let mut writer = Vec::with_capacity(128);
        self.to_map_writer(&mut writer, reads, mappings).unwrap();
        String::from_utf8(writer).unwrap()
    }
    ///
    /// create MAP file
    ///
    pub fn to_map_writer<W: std::io::Write, S: Seq>(
        &self,
        mut writer: W,
        reads: &ReadCollection<S>,
        mappings: &Mappings,
    ) -> std::io::Result<()> {
        // header
        writeln!(writer, "# read\tpos\tbase\tnodes_and_probs")?;

        // body
        for (i, read) in reads.iter().enumerate() {
            for (j, &base) in read.as_ref().into_iter().enumerate() {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}",
                    i,
                    j,
                    base as char,
                    mappings[i].to_nodes_string(j)
                )?;
            }
        }
        Ok(())
    }
    ///
    ///
    ///
    pub fn from_map_reader<R: std::io::BufRead, S: Seq>(
        &self,
        reader: R,
        reads: &ReadCollection<S>,
    ) -> Mappings {
        let mut ret: Vec<Vec<Vec<(NodeIndex, Prob)>>> = vec![];
        for (_, read) in reads.into_iter().enumerate() {
            ret.push(vec![vec![]; read.as_ref().len()]);
        }
        for line in reader.lines() {
            let text = line.unwrap();
            let first_char = text.chars().nth(0).unwrap();
            if first_char != '#' {
                let mut iter = text.split_whitespace();

                let i: usize = iter.next().unwrap().parse().unwrap();
                let j: usize = iter.next().unwrap().parse().unwrap();
                iter.next().unwrap(); // base
                let p: Vec<(NodeIndex, Prob)> = iter
                    .next()
                    .unwrap()
                    .split(',')
                    .map(|s| {
                        let mut sp = s.split(':');
                        let v: usize = sp.next().unwrap().parse().unwrap();
                        let p: f64 = sp.next().unwrap().parse().unwrap();
                        (NodeIndex::new(v), Prob::from_prob(p / 100.0))
                    })
                    .collect();
                ret[i][j] = p;
            }
        }

        let mut ms = vec![];
        for (i, _read) in reads.into_iter().enumerate() {
            let m = Mapping::from_nodes_and_probs(&ret[i]);
            ms.push(m);
        }
        Mappings::new(ms)
    }
    ///
    /// parse MAP string with [`MultiDbg::from_map_reader`]
    ///
    pub fn from_map_str<S: Seq>(&self, s: &str, reads: &ReadCollection<S>) -> Mappings {
        self.from_map_reader(s.as_bytes(), reads)
    }
    ///
    /// parse MAP file with [`MultiDbg::from_map_reader`]
    ///
    pub fn from_map_file<P: AsRef<std::path::Path>, S: Seq>(
        &self,
        path: P,
        reads: &ReadCollection<S>,
    ) -> Mappings {
        let file = std::fs::File::open(path).unwrap();
        let reader = std::io::BufReader::new(file);
        self.from_map_reader(reader, reads)
    }
}

///
/// GFA format
///
impl MultiDbg {
    /// GFA format
    ///
    /// ```text
    /// # comment
    ///
    /// # segment
    /// #  name  sequence   copy_number  length
    /// S  0     ATCGATTCG  DP:i:1       LN:i:9
    /// S  1     ATCGATTCG  DP:i:1       LN:i:9
    ///
    /// # link
    /// #  source  orient  target  orient  optional  node_id
    /// L  0       +       1       +       *         ID:Z:0
    /// ```
    ///
    /// * segment for each edge
    /// * link from in_edge to out_edge of each node
    ///
    pub fn to_gfa_writer_with<W, FL, FC>(
        &self,
        mut writer: W,
        label: FL,
        color: FC,
    ) -> std::io::Result<()>
    where
        W: std::io::Write,
        FL: Fn(EdgeIndex) -> String,
        FC: Fn(EdgeIndex) -> (u8, u8, u8),
    {
        for (edge, s, t, weight) in self.edges_compact() {
            let seq = &self.seq_compact(edge);
            let (r, g, b) = color(edge);
            let color = format!("#{:02x}{:02x}{:02x}", r, g, b);
            writeln!(
                writer,
                "S\t{}\t{}\tDP:f:{}\tLN:i:{}\tLB:Z:{}\tCL:Z:{}",
                edge.index(),
                sequence_to_string(&seq),
                self.copy_num_of_edge_in_compact(edge),
                seq.len(),
                label(edge),
                color,
            )?
        }
        let terminal = self.terminal_node_compact();
        for (node, weight) in self.nodes_compact() {
            if terminal.is_some() && node != terminal.unwrap() {
                for (in_edge, _, _) in self.parents_compact(node) {
                    for (out_edge, _, _) in self.childs_compact(node) {
                        writeln!(
                            writer,
                            "L\t{}\t+\t{}\t+\t*\tID:Z:{}",
                            in_edge.index(),
                            out_edge.index(),
                            node.index(),
                        )?
                    }
                }
            }
        }
        Ok(())
    }
    /// Default GFA format with sequence in label
    ///
    pub fn to_gfa_writer<W: std::io::Write>(&self, writer: W) -> std::io::Result<()> {
        self.to_gfa_writer_with(
            writer,
            |e| sequence_to_string(&self.seq_compact(e)).to_string(),
            |e| (0, 255, 0),
        )
    }
    ///
    /// create GFA string with `to_gfa_writer`
    ///
    pub fn to_gfa_string(&self) -> String {
        let mut writer = Vec::with_capacity(128);
        self.to_gfa_writer(&mut writer).unwrap();
        String::from_utf8(writer).unwrap()
    }
    ///
    /// create GFA file with `to_gfa_writer`
    ///
    pub fn to_gfa_file<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        let mut file = std::fs::File::create(path).unwrap();
        self.to_gfa_writer(&mut file)
    }
    ///
    /// GFA-post format with posterior information in LB
    ///
    pub fn to_gfa_post_writer<W: std::io::Write>(
        &self,
        writer: W,
        posterior: &Posterior,
        copy_nums_true: Option<&CopyNums>,
    ) -> std::io::Result<()> {
        self.to_gfa_writer_with(
            writer,
            |e| {
                let p_edge = posterior.p_edge(e);
                match copy_nums_true {
                    Some(copy_nums_true) => {
                        format!(
                            "{:.2}x,{}x({})",
                            p_edge.mean(),
                            copy_nums_true[e],
                            p_edge.to_short_string()
                        )
                    }
                    None => {
                        format!("{:.2}x,?x({})", p_edge.mean(), p_edge.to_short_string())
                    }
                }
            },
            |e| match copy_nums_true {
                Some(copy_nums_true) => {
                    let copy_num = posterior.p_edge(e).mean();
                    let copy_num_true = copy_nums_true[e] as f64;
                    let max = 200_u8;
                    let half = (max / 2) as f64;

                    if copy_num > copy_num_true {
                        // over-represented: red
                        let r = (((copy_num - copy_num_true) * half) as u8).clamp(0, max);
                        (max, max - r, max - r)
                    } else {
                        // under-represented: blue
                        let b = (((copy_num_true - copy_num) * half) as u8).clamp(0, max);
                        (max - b, max - b, max)
                    }
                }
                None => (0, 0, 0),
            },
        )
    }
    ///
    /// create GFA-post file with posterior information
    ///
    pub fn to_gfa_post_file<P: AsRef<std::path::Path>>(
        &self,
        path: P,
        posterior: &Posterior,
        copy_nums_true: Option<&CopyNums>,
    ) -> std::io::Result<()> {
        let mut file = std::fs::File::create(path).unwrap();
        self.to_gfa_post_writer(&mut file, posterior, copy_nums_true)
    }
}

impl std::fmt::Display for MultiDbg {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let seqs = self
            .to_styled_seqs()
            .iter()
            .map(|seq| seq.to_string())
            .join(",");
        write!(f, "{},{}", self.k(), seqs)
    }
}

impl std::fmt::Display for MultiCompactEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.edges_in_full
                .iter()
                .map(|e| format!("e{}", e.index()))
                .join(",")
        )
    }
}

impl std::fmt::Display for MultiFullEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}({}x)", self.base as char, self.copy_num)
    }
}

impl std::fmt::Display for MultiCompactNode {
    fn fmt(&self, _: &mut std::fmt::Formatter) -> std::fmt::Result {
        Ok(())
    }
}

impl std::fmt::Display for MultiFullNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "is_terminal={}", self.is_terminal)
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::super::toy;
    use super::*;
    use crate::hmmv2::params::PHMMParams;

    fn assert_dbg_dumpload_is_correct(dbg: MultiDbg) {
        dbg.show_graph_with_kmer();

        let s = dbg.to_dbg_string();
        println!("{}", s);

        dbg.to_dbg_file("hoge.dbg");

        let dbg_1 = MultiDbg::from_dbg_str(&s);
        dbg_1.show_graph_with_kmer();
        let s_1 = dbg.to_dbg_string();

        assert!(dbg.is_equivalent(&dbg_1));
        assert!(s == s_1);
    }
    #[test]
    fn dumpload() {
        assert_dbg_dumpload_is_correct(toy::circular());
        assert_dbg_dumpload_is_correct(toy::linear());
        assert_dbg_dumpload_is_correct(toy::intersection());
        assert_dbg_dumpload_is_correct(toy::selfloop());
        assert_dbg_dumpload_is_correct(toy::repeat());
    }
    #[test]
    fn gfa() {
        let dbg = toy::repeat();
        dbg.show_graph_with_kmer();
        let s = dbg.to_gfa_string();
        println!("{}", s);

        dbg.to_gfa_file("repeat.gfa");
    }
    #[test]
    fn dbg_gz_compressed() {
        let dbg_a = toy::repeat();
        dbg_a.show_graph_with_kmer();
        dbg_a.to_dbg_file("repeat.gfa.gz");

        let dbg_b = MultiDbg::from_dbg_file("repeat.gfa.gz");
        dbg_b.show_graph_with_kmer();
        // assert!(dbg_b.is_equal(&dbg_a));
    }
    #[test]
    fn map() {
        let dbg = toy::intersection();
        dbg.show_graph_with_kmer();
        let reads = ReadCollection {
            reads: vec![b"GATCC".to_vec(), b"TATCA".to_vec()],
        };
        let param = PHMMParams::default();
        let mappings = dbg.generate_mappings(param, &reads, None);
        let s = dbg.to_map_string(&reads, &mappings);
        println!("{}", s);

        let mappings2 = dbg.from_map_str(&s, &reads);
        println!("{}", mappings2[0]);
        println!("{}", mappings2[1]);
        for (i, read) in reads.into_iter().enumerate() {
            for j in 0..read.len() {
                assert_eq!(mappings[i].nodes(j), mappings2[i].nodes(j));
            }
        }
    }
}
