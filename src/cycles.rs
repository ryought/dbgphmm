use crate::dbg::{DbgHash, DBG};
use crate::kmer::kmer::Kmer;
use fnv::FnvHashMap as HashMap;
use fnv::FnvHashSet as HashSet;
use log::{debug, info, warn};
use std::collections::VecDeque;

pub struct DbgTree {
    root: Kmer,
    parent_on_tree: HashMap<Kmer, Kmer>,
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
        let mut parent_on_tree: HashMap<Kmer, Kmer> = HashMap::default();
        let mut depth: HashMap<Kmer, usize> = HashMap::default();

        // do dfs
        let mut is_visited: HashSet<Kmer> = HashSet::default();
        let mut q: VecDeque<(Kmer, usize)> = VecDeque::new();

        // start from root
        is_visited.insert(root.prefix());
        is_visited.insert(root.suffix());
        depth.insert(root.clone(), 0);
        q.push_back((root.clone(), 0));

        while !q.is_empty() {
            let (kmer, d) = q.pop_front().unwrap();
            debug!("DFS: now: {}@{}", kmer, d);

            for child in dbg.childs(&kmer).into_iter() {
                let suffix = child.suffix();
                if !is_visited.contains(&suffix) {
                    // kmer -> child move
                    debug!("DFS: hop_child: {} -> {}", kmer, child);
                    is_visited.insert(suffix);
                    parent_on_tree.insert(child.clone(), kmer.clone());
                    depth.insert(child.clone(), d + 1);
                    q.push_back((child, d + 1));
                }
            }

            for parent in dbg.parents(&kmer).into_iter() {
                let prefix = parent.prefix();
                if !is_visited.contains(&prefix) {
                    // kmer -> parent move
                    debug!("DFS: hop_parent: {} -> {}", kmer, parent);
                    is_visited.insert(prefix);
                    parent_on_tree.insert(parent.clone(), kmer.clone());
                    depth.insert(parent.clone(), d + 1);
                    q.push_back((parent, d + 1));
                }
            }
        }

        let mut loop_edges: Vec<Kmer> = Vec::new();
        for kmer in dbg.kmers().into_iter() {
            if !depth.contains_key(&kmer) {
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
    pub fn cycle_components(&self, kmer: &Kmer) -> Vec<Kmer> {
        let mut lefts: Vec<Kmer> = Vec::new();
        let mut rights: Vec<Kmer> = Vec::new();

        // find joints in the tree
        let (mut left, left_depth) = kmer
            .childs()
            .into_iter()
            .find_map(|child| match self.depth.get(&child) {
                Some(&depth) => Some((child, depth)),
                None => None,
            })
            .unwrap();
        let (mut right, right_depth) = kmer
            .parents()
            .into_iter()
            .find_map(|parent| match self.depth.get(&parent) {
                Some(&depth) => Some((parent, depth)),
                None => None,
            })
            .unwrap();
        debug!("{}@{} -- {}@{}", left, left_depth, right, right_depth);

        // move to same depth
        let d = left_depth.min(right_depth);
        debug!("depth: {}", d);

        for i in d..right_depth {
            let new_right = self.parent_on_tree.get(&right).unwrap().clone();
            rights.push(right);
            right = new_right;
        }
        for i in d..left_depth {
            let new_left = self.parent_on_tree.get(&left).unwrap().clone();
            lefts.push(left);
            left = new_left;
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
                rights.push(right);
                break;
            } else {
                let new_right = self.parent_on_tree.get(&right).unwrap().clone();
                rights.push(right);
                right = new_right;
                let new_left = self.parent_on_tree.get(&left).unwrap().clone();
                lefts.push(left);
                left = new_left;
            }
        }

        let path: Vec<Kmer> = std::iter::once(kmer.clone())
            .chain(rights.into_iter())
            .chain(lefts.into_iter().rev())
            .collect();
        path
    }
}
