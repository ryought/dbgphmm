//!
//! Cycle enumeration using [Johnson1975]() algorithm
//!
//! both for directed and undirected graphs
//!
//! * Tarjan1972
//! Enumeration of the Elementary Circuits of a Directed Graph
//! https://ecommons.cornell.edu/handle/1813/5941
//!
//! * Johnson1975
//! Finding all the elementary circuits of a directed graph
//!
use fnv::FnvHashSet as HashSet;
use petgraph::graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph};
use petgraph::Direction;

///
/// Enumerate all simple cycles (elementary circuits; a node-simple path starts/ends from the same node)
///
/// Assumptions
/// * graph is strongly-connected.
///     there are path of both directions (v->w and w->v) between any two nodes v and w).
///
/// TODO
/// * if graph has self-loops, output them as simple cycles.
///
pub fn simple_cycles<N, E>(graph: &DiGraph<N, E>) -> Vec<Vec<NodeIndex>> {
    // let mut cycles = vec![];
    let n = graph.node_count();
    let ix = |node: NodeIndex| node.index();
    let mut ret: Vec<Vec<NodeIndex>> = Vec::new();
    for start_node in graph.node_indices() {
        // search for circuit starting from node
        //
        let mut blocked = vec![false; n];
        blocked[ix(start_node)] = true;
        // B-list
        let mut b: Vec<HashSet<NodeIndex>> = vec![HashSet::default(); n];
        // abbr of ix() node.index()
        // recursively unblocks nodes connected by B-list
        let unblock =
            |node: NodeIndex, blocked: &mut Vec<bool>, b: &mut Vec<HashSet<NodeIndex>>| {
                let mut nodes = vec![node];
                while let Some(node) = nodes.pop() {
                    blocked[ix(node)] = false;
                    for &bnode in b[ix(node)].iter() {
                        nodes.push(bnode);
                    }
                    b[ix(node)].clear();
                }
            };
        let neighbors = |node: NodeIndex| {
            let ret: Vec<_> = graph
                .neighbors_directed(node, Direction::Outgoing)
                .filter(|next_node| next_node.index() >= start_node.index())
                .collect();
            ret
        };
        let mut path: Vec<NodeIndex> = vec![start_node];
        let mut stack: Vec<(NodeIndex, Vec<NodeIndex>)> = vec![(start_node, neighbors(start_node))];
        let mut closed: HashSet<NodeIndex> = HashSet::default();

        // CIRCUIT(node) routine
        while let Some((node, next_nodes)) = stack.last_mut() {
            // for w in A_k(v)
            if !next_nodes.is_empty() {
                let next_node = next_nodes.pop().unwrap();
                // there are neighbor nodes
                // (1) back to the start node
                if next_node == start_node {
                    // cycle found!
                    ret.push(path.clone());
                    for &node_in_path in path.iter() {
                        closed.insert(node_in_path);
                    }
                }
                // (2) visit a unblocked neighboring node
                if !blocked[ix(next_node)] {
                    path.push(next_node);
                    blocked[ix(next_node)] = true;
                    stack.push((next_node, neighbors(next_node)));
                    closed.remove(&next_node);
                    continue;
                }
            }
            // if f then unblock else put v in b(w)
            if next_nodes.is_empty() {
                // no neighbors
                if closed.contains(node) {
                    // f=true
                    unblock(*node, &mut blocked, &mut b);
                } else {
                    // f=false
                    for neighbor in neighbors(*node) {
                        b[ix(neighbor)].insert(*node);
                    }
                }
                stack.pop();
                path.pop();
            }
        }
    }
    ret
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;

    #[test]
    fn simple_cycles_test_1() {
        let g: DiGraph<(), ()> =
            DiGraph::from_edges(&[(0, 0), (0, 1), (0, 2), (1, 2), (2, 0), (2, 1), (2, 2)]);
        let cycles = simple_cycles(&g);
        println!("cycles={:?}", cycles);
        assert_eq!(
            cycles,
            vec![
                vec![ni(0)],
                vec![ni(0), ni(1), ni(2)],
                vec![ni(0), ni(2)],
                vec![ni(1), ni(2)],
                vec![ni(2)],
            ]
        );
    }
}
