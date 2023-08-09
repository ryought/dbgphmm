//!
//! Posterior output formats (INSPECT format)
//!
use super::{Posterior, PosteriorSample, Score};
use crate::multi_dbg::{CopyNums, MultiDbg};
use crate::prob::Prob;
use itertools::Itertools;

///
/// dump and load functions
///
/// ```text
/// Z   -19281.0228
/// C   -192919.0    [1,2,1,1,1,2,1,0]   likelihood=0.00 e1+e2-e3+,e1+e2-
/// C   -191882.0    [1,2,0,0,1,2,2,1]   likelihood=0.01
/// ```
///
impl Posterior {
    ///
    ///
    ///
    pub fn to_writer<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        writeln!(writer, "# {}", env!("GIT_HASH"))?;
        writeln!(writer, "Z\t{}", self.p.to_log_value())?;
        for sample in self
            .samples
            .iter()
            .sorted_by_key(|sample| sample.score.p())
            .rev()
        {
            writeln!(
                writer,
                "C\t{}\t{}\t{}\t{}",
                sample.score.p().to_log_value(),
                sample.copy_nums,
                sample.score,
                sample.to_infos_string(),
            )?
        }
        Ok(())
    }
    ///
    /// create string with `to_gfa_writer`
    ///
    pub fn to_string(&self) -> String {
        let mut writer = Vec::with_capacity(128);
        self.to_writer(&mut writer).unwrap();
        String::from_utf8(writer).unwrap()
    }
    ///
    /// create file with `to_gfa_writer`
    ///
    pub fn to_file<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        let mut file = std::fs::File::create(path).unwrap();
        self.to_writer(&mut file)
    }
    ///
    ///
    ///
    pub fn from_reader<R: std::io::BufRead>(reader: R) -> Self {
        let mut samples = Vec::new();
        let mut p = Prob::zero();

        for line in reader.lines() {
            let text = line.unwrap();
            let first_char = text.chars().nth(0).unwrap();
            match first_char {
                'C' => {
                    let mut iter = text.split_whitespace();
                    iter.next().unwrap(); // 'C'
                    iter.next().unwrap(); // value
                    let copy_nums: CopyNums = iter.next().unwrap().parse().unwrap();
                    let score: Score = iter.next().unwrap().parse().unwrap();
                    let infos = PosteriorSample::from_infos_str(iter.next().unwrap()).unwrap();
                    p += score.p();
                    samples.push(PosteriorSample {
                        copy_nums,
                        score,
                        infos,
                    });
                }
                _ => {} // ignore
            }
        }

        Posterior { samples, p }
    }
    ///
    ///
    pub fn from_str(s: &str) -> Self {
        Self::from_reader(s.as_bytes())
    }
    ///
    ///
    pub fn from_file<P: AsRef<std::path::Path>>(path: P) -> Self {
        let file = std::fs::File::open(path).unwrap();
        let reader = std::io::BufReader::new(file);
        Self::from_reader(reader)
    }
}

fn format_option_copy_num<T: std::fmt::Display>(o: Option<T>) -> String {
    match o {
        Some(i) => format!("{}", i),
        None => format!("?"),
    }
}

///
/// benchmark functions for when true genome is available
///
impl MultiDbg {
    ///
    /// Everytime
    /// * posterior probability (normalized)
    /// * likelihood (log)
    /// * prior (log)
    /// * genome size
    ///
    /// Only if genome is known
    /// * diff of copynums from true
    ///
    pub fn to_inspect_writer<W: std::io::Write>(
        &self,
        mut writer: W,
        posterior: &Posterior,
        copy_nums_true: Option<&CopyNums>,
    ) -> std::io::Result<()> {
        // for each copy nums
        // TODO clean up using key value format
        writeln!(writer, "# {}", env!("GIT_HASH"))?;
        writeln!(
            writer,
            "{}\tG\tn_edges_full\t{}",
            self.k(),
            self.n_edges_full()
        )?;
        writeln!(
            writer,
            "{}\tG\tn_edges_compact\t{}",
            self.k(),
            self.n_edges_compact()
        )?;
        writeln!(
            writer,
            "{}\tG\tn_nodes_full\t{}",
            self.k(),
            self.n_nodes_full()
        )?;
        writeln!(
            writer,
            "{}\tG\tn_nodes_compact\t{}",
            self.k(),
            self.n_nodes_compact()
        )?;
        writeln!(
            writer,
            "{}\tG\tn_emittable_edges\t{}",
            self.k(),
            self.n_emittable_edges()
        )?;
        writeln!(
            writer,
            "{}\tG\tdegree_stats\t{:?}",
            self.k(),
            self.degree_stats(),
        )?;
        for (i, sample) in posterior
            .samples
            .iter()
            .sorted_by_key(|sample| sample.score.p())
            .rev()
            .enumerate()
        {
            let score = &sample.score;
            let copy_nums = &sample.copy_nums;
            let json = serde_json::to_string(score)?;
            writeln!(
                writer,
                "{}\tC\t{}\t{:.10}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.k(),
                i,
                (score.p() / posterior.p()).to_value(),
                score.likelihood.to_log_value(),
                score.prior.to_log_value(),
                score.n_euler_circuits,
                score.genome_size,
                format_option_copy_num(
                    copy_nums_true.map(|copy_nums_true| copy_nums_true.diff(copy_nums))
                ),
                sample.to_infos_string(),
                copy_nums,
                json,
            )?
        }

        // for each edges
        for edge in self.graph_compact().edge_indices() {
            let p_edge = posterior.p_edge(edge);
            let copy_num_true = copy_nums_true.map(|copy_num_true| copy_num_true[edge]);
            writeln!(
                writer,
                "{}\tE\te{}\t{}\t{:.5}\t{:.5}\t{:.5}\t{}",
                self.k(),
                edge.index(),
                format_option_copy_num(copy_num_true),
                p_edge.mean(),
                format_option_copy_num(
                    copy_num_true.map(|copy_num_true| p_edge.p_x(copy_num_true).to_value())
                ),
                p_edge.p_x(0).to_value(),
                p_edge.to_short_string(),
            )?
        }

        Ok(())
    }
    ///
    ///
    pub fn to_inspect_string(
        &self,
        posterior: &Posterior,
        copy_nums_true: Option<&CopyNums>,
    ) -> String {
        let mut writer = Vec::with_capacity(128);
        self.to_inspect_writer(&mut writer, posterior, copy_nums_true)
            .unwrap();
        String::from_utf8(writer).unwrap()
    }
    ///
    ///
    pub fn to_inspect_file<P: AsRef<std::path::Path>>(
        &self,
        path: P,
        posterior: &Posterior,
        copy_nums_true: Option<&CopyNums>,
    ) -> std::io::Result<()> {
        let mut file = std::fs::File::create(path).unwrap();
        self.to_inspect_writer(&mut file, posterior, copy_nums_true)
    }
    ///
    /// Parse INSPECT file
    ///
    pub fn from_inspect_reader<R: std::io::BufRead>(&self, reader: R) -> Option<Posterior> {
        let mut posterior = Posterior::new();

        for line in reader.lines() {
            let text = line.unwrap();
            let first_char = text.chars().nth(0)?;
            if first_char == '#' {
                continue;
            }
            let columns: Vec<&str> = text.split_whitespace().collect();

            let k: usize = columns.get(0)?.parse().ok()?;
            assert_eq!(k, self.k());

            let kind = columns.get(1)?.chars().nth(0)?;
            match kind {
                'C' => {
                    let infos = PosteriorSample::from_infos_str(columns.get(9)?)?;
                    let copy_nums = columns.get(10)?.parse().ok()?;
                    let score: Score = serde_json::from_str(columns.get(11)?).ok()?;
                    posterior.add(PosteriorSample {
                        copy_nums,
                        score,
                        infos,
                    });
                }
                _ => {
                    // TODO
                    // parse true copy_nums?
                }
            }
        }

        Some(posterior)
    }
    ///
    ///
    pub fn from_inspect_str(&self, s: &str) -> Option<Posterior> {
        self.from_inspect_reader(s.as_bytes())
    }
    ///
    ///
    pub fn from_inspect_file<P: AsRef<std::path::Path>>(&self, path: P) -> Option<Posterior> {
        let file = std::fs::File::open(path.as_ref()).unwrap();
        let reader = std::io::BufReader::new(file);
        self.from_inspect_reader(reader)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ei;
    use crate::multi_dbg::neighbors::{UpdateInfo, UpdateMethod};
    use crate::multi_dbg::toy;
    use crate::prob::p;
    use rustflow::min_flow::residue::ResidueDirection;

    fn example_posterior() -> (MultiDbg, Posterior) {
        let mut dbg = toy::intersection();
        let mut post = Posterior::new();
        post.add(PosteriorSample {
            copy_nums: vec![1, 1, 2, 2].into(),
            score: Score {
                likelihood: p(0.6),
                prior: p(0.3),
                time_likelihood: 10,
                genome_size: 101,
                n_euler_circuits: 10.0,
                time_euler: 12,
            },
            infos: vec![
                UpdateInfo::new(
                    vec![vec![
                        (ei(0), ResidueDirection::Up),
                        (ei(121), ResidueDirection::Down),
                    ]],
                    UpdateMethod::Manual,
                ),
                UpdateInfo::new(
                    vec![vec![
                        (ei(3), ResidueDirection::Down),
                        (ei(1), ResidueDirection::Up),
                    ]],
                    UpdateMethod::Manual,
                ),
            ],
        });
        post.add(PosteriorSample {
            copy_nums: vec![1, 1, 1, 2].into(),
            score: Score {
                likelihood: p(0.003),
                prior: p(0.2),
                time_likelihood: 11,
                genome_size: 99,
                n_euler_circuits: 10.0,
                time_euler: 12,
            },
            infos: Vec::new(),
        });
        (dbg, post)
    }

    #[test]
    fn posterior_dump_load() {
        let (_, post) = example_posterior();
        let s = post.to_string();
        println!("{}", s);

        let post_loaded = Posterior::from_str(&s);
        assert_eq!(post_loaded.samples, post.samples);
        assert_eq!(post_loaded.p, post.p);
    }

    #[test]
    fn parse_inspect() {
        // parsing Vec<UpdateCycle> = Vec<Vec<(EdgeIndex, ResidueDirection)>>

        // case 1:
        // containing multiple updates
        let infos = vec![
            UpdateInfo::new(
                vec![vec![
                    (ei(5), ResidueDirection::Up),
                    (ei(6), ResidueDirection::Up),
                    (ei(1), ResidueDirection::Down),
                ]],
                UpdateMethod::Manual,
            ),
            UpdateInfo::new(
                vec![vec![
                    (ei(10), ResidueDirection::Up),
                    (ei(1), ResidueDirection::Up),
                ]],
                UpdateMethod::Rescue {
                    index: 0,
                    n_kmers: 50,
                    length: 111,
                    freq: 5.001,
                    non_zero: false,
                },
            ),
            UpdateInfo::new(
                vec![vec![(ei(5), ResidueDirection::Down)]],
                UpdateMethod::Short,
            ),
        ];
        let s = "[X(e5+e6+e1-),R(e10+e1+|0|50|111|5.001|Z),S(e5-)]";
        let infos2 = PosteriorSample::from_infos_str(s);
        assert_eq!(PosteriorSample::to_infos_string_internal(&infos), s);
        assert_eq!(infos2, Some(infos));

        // case 2:
        // containing no cycle (or update)
        let infos = PosteriorSample::from_infos_str("[]");
        assert_eq!(infos, Some(vec![]));
        assert_eq!(PosteriorSample::to_infos_string_internal(&[]), "[]");

        // full INSPECT file
        let (dbg, post) = example_posterior();
        let s = dbg.to_inspect_string(&post, None);
        println!("{}", s);
        let post2 = dbg.from_inspect_str(&s).unwrap();
        println!("{:?}", post2);
        assert_eq!(post.samples, post2.samples);
    }
}
