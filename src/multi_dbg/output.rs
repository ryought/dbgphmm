//!
//! Output (serialization and deserialization) functions of MultiDbg
//!
use super::{MultiCompactEdge, MultiCompactNode, MultiDbg, MultiFullEdge, MultiFullNode};
use crate::common::{sequence_to_string, CopyNum, NULL_BASE};
use itertools::Itertools;
use petgraph::graph::{DefaultIx, DiGraph, EdgeIndex, NodeIndex};

///
/// serialize/deserialize and debug print methods
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
        let mut file = std::fs::File::create(path).unwrap();
        self.to_dbg_writer(&mut file)
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
                _ => panic!("invalid DBG format"),
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
        for _ in nodes.iter() {
            compact.add_node(MultiCompactNode::new());
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
        let file = std::fs::File::open(path).unwrap();
        let reader = std::io::BufReader::new(file);
        Self::from_dbg_reader(reader)
    }
    /// GFA format
    ///
    /// ```text
    /// # comment
    ///
    /// # segment
    /// #  name  sequence   copy_number
    /// S  0     ATCGATTCG  CN:i:1
    /// S  1     ATCGATTCG  CN:i:1
    ///
    /// # link
    /// #  source  orient  target  orient  optional  node_id
    /// L  0       +       1       +       *         ID:Z:0
    /// ```
    ///
    /// * segment for each edge
    /// * link from in_edge to out_edge of each node
    ///
    pub fn to_gfa_writer<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        for (edge, s, t, weight) in self.edges_compact() {
            let seq = &self.seq_compact(edge);
            writeln!(
                writer,
                "S\t{}\t{}\tDP:f:{}\tLB:Z:{}\tLN:i:{}",
                edge.index(),
                sequence_to_string(&seq),
                self.copy_num_of_edge_in_compact(edge),
                sequence_to_string(&seq),
                seq.len(),
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
}
