use crate::dbg::{DbgHash, DBG};
use crate::kmer::kmer::Kmer;
use fnv::FnvHashMap as HashMap;
use fnv::FnvHashSet as HashSet;
use log::{debug, info, warn};
use std::collections::VecDeque;
use std::fmt::Write as FmtWrite;

#[derive(Copy, Clone)]
enum ParentType {
    /// parent XYYYY -> child YYYYZ
    Succ,
    /// parent YYYYX -> child ZYYYY
    Pred,
}

/// Cycle is a sequence of k-mers
/// For example, direction of ATCG -> TCGT -> CGTA is all forward.
/// but ATCG -> CATC is reverse.
/// i.e.
/// Forward <=> XYYY -> YYYZ
/// Reverse <=> YYYX -> ZYYY
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum CycleDirection {
    Forward,
    Reverse,
}

pub struct DbgTree {
    root: Kmer,
    parent_on_tree: HashMap<Kmer, (Kmer, ParentType)>,
    depth: HashMap<Kmer, usize>,
    loop_edges: Vec<Kmer>,
}

///
/// Store cycles and spanning tree.
/// fast cycle lookup on DBG.
///
impl DbgTree {
    pub fn new<T: DBG>(dbg: &T, root: &Kmer) -> DbgTree {
        // struct elements
        let mut parent_on_tree: HashMap<Kmer, (Kmer, ParentType)> = HashMap::default(); // parent_on_tree[node] = (node, succ|pred)
        let mut depth: HashMap<Kmer, usize> = HashMap::default(); // depth of node

        // do dfs
        let mut is_visited: HashSet<Kmer> = HashSet::default(); // node is visited or not
        let mut is_used_edge: HashSet<Kmer> = HashSet::default(); // edge is used or not
        let mut q: VecDeque<(Kmer, usize)> = VecDeque::new(); // bfs queue (node, depth)

        // start from root
        is_visited.insert(root.clone());
        depth.insert(root.clone(), 0);
        q.push_back((root.clone(), 0));

        while !q.is_empty() {
            let (km1mer, d) = q.pop_front().unwrap();
            debug!("DFS: now: {}@{}", km1mer, d);

            // km1mer XYYY -> succ YYYZ
            for succ in dbg.succs(&km1mer).into_iter() {
                let suffix = succ.suffix();
                if !is_visited.contains(&suffix) {
                    debug!("DFS: hop_succ: {} -> {}", km1mer, suffix);
                    // suffix=YYZ <- km1mer=YYY by moving edge=YYYZ
                    parent_on_tree.insert(suffix.clone(), (km1mer.clone(), ParentType::Succ));
                    is_visited.insert(suffix.clone());
                    is_used_edge.insert(succ);
                    depth.insert(suffix.clone(), d + 1);
                    q.push_back((suffix, d + 1));
                }
            }

            // km1mer YYYX -> preds ZYYY
            for pred in dbg.preds(&km1mer).into_iter() {
                let prefix = pred.prefix();
                if !is_visited.contains(&prefix) {
                    debug!("DFS: hop_pred: {} -> {}", km1mer, prefix);
                    // prefix=ZYY <- km1mer=YYY by moving edge=ZYYY
                    parent_on_tree.insert(prefix.clone(), (km1mer.clone(), ParentType::Pred));
                    is_visited.insert(prefix.clone());
                    is_used_edge.insert(pred);
                    depth.insert(prefix.clone(), d + 1);
                    q.push_back((prefix, d + 1));
                }
            }
        }

        // check if the tree is spanning whole graph?
        for kmer in dbg.kmers().iter() {
            let suffix = kmer.suffix();
            let prefix = kmer.prefix();
            if !is_visited.contains(&suffix) {
                warn!("not spanning: {}", suffix);
            }
            if !is_visited.contains(&prefix) {
                warn!("not spanning: {}", prefix);
            }
        }

        let mut loop_edges: Vec<Kmer> = Vec::new();
        for kmer in dbg.kmers().into_iter() {
            if !is_used_edge.contains(&kmer) {
                loop_edges.push(kmer);
            }
        }

        DbgTree {
            root: root.clone(),
            parent_on_tree,
            depth,
            loop_edges,
        }
    }
    pub fn root(&self) -> &Kmer {
        &self.root
    }
    pub fn cycle_keys(&self) -> &[Kmer] {
        &self.loop_edges[..]
    }
    fn go_up_tree(&self, km1mer: &Kmer) -> (Kmer, Kmer, ParentType) {
        let (parent_km1mer, parent_type) = self.parent_on_tree.get(km1mer).unwrap().clone();
        let edge_kmer = match parent_type {
            ParentType::Succ => km1mer.extend_first(parent_km1mer.first()),
            ParentType::Pred => km1mer.extend_last(parent_km1mer.last()),
        };
        (parent_km1mer.clone(), edge_kmer, parent_type)
    }
    pub fn cycle_components(&self, kmer: &Kmer) -> Vec<(Kmer, CycleDirection)> {
        let mut lefts: Vec<(Kmer, CycleDirection)> = Vec::new(); // edges
        let mut rights: Vec<(Kmer, CycleDirection)> = Vec::new(); // edges

        // find joints in the tree
        //  <------   left
        //   XXXXXXX
        //    ------> right
        let mut left = kmer.prefix();
        let &left_depth = self.depth.get(&left).unwrap_or_else(|| panic!("left!"));
        let mut right = kmer.suffix();
        let &right_depth = self.depth.get(&right).unwrap_or_else(|| panic!("right!"));
        debug!("{}@{} -- {}@{}", left, left_depth, right, right_depth);

        // move to same depth
        let d = left_depth.min(right_depth);
        debug!("depth: {}", d);

        for _ in d..right_depth {
            debug!("right: {}", d);
            let (r, edge, p_type) = self.go_up_tree(&right);
            // if going right (child -> parent), Pred is Forward
            let direction = match p_type {
                ParentType::Pred => CycleDirection::Forward,
                ParentType::Succ => CycleDirection::Reverse,
            };
            rights.push((edge, direction));
            right = r;
        }
        for _ in d..left_depth {
            debug!("left: {}", d);
            let (l, edge, p_type) = self.go_up_tree(&left);
            // if going left (parent -> child), Succ is Forward
            let direction = match p_type {
                ParentType::Succ => CycleDirection::Forward,
                ParentType::Pred => CycleDirection::Reverse,
            };
            lefts.push((edge, direction));
            left = l;
        }

        // find LCA
        for i in 0..=d {
            debug!(
                "same: {}@{} {}@{}",
                left,
                self.depth.get(&left).unwrap(),
                right,
                self.depth.get(&right).unwrap()
            );
            if right == left || i == d {
                break;
            } else {
                let (r, edge, p_type) = self.go_up_tree(&right);
                // if going right (child -> parent), Pred is Forward
                let direction = match p_type {
                    ParentType::Pred => CycleDirection::Forward,
                    ParentType::Succ => CycleDirection::Reverse,
                };
                rights.push((edge, direction));
                right = r;
                let (l, edge, p_type) = self.go_up_tree(&left);
                // if going left (parent -> child), Succ is Forward
                let direction = match p_type {
                    ParentType::Succ => CycleDirection::Forward,
                    ParentType::Pred => CycleDirection::Reverse,
                };
                lefts.push((edge, direction));
                left = l;
            }
        }

        let path: Vec<(Kmer, CycleDirection)> =
            std::iter::once((kmer.clone(), CycleDirection::Forward))
                .chain(rights.into_iter())
                .chain(lefts.into_iter().rev())
                .collect();
        path
    }
    pub fn as_stats(&self) -> String {
        let mut s = String::new();
        writeln!(&mut s, "#LoopEdge: {}", self.loop_edges.len());
        // TODO cycle length histogram?
        for (i, e) in self.loop_edges.iter().enumerate() {
            writeln!(
                &mut s,
                "cycle #{} {} {}",
                i,
                e,
                self.cycle_components(e).len()
            );
            /*
            for p in s.cycle_components(e).iter() {
                info!("{}", p);
            }
            */
        }
        s
    }
}
