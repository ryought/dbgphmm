use fnv::FnvHashMap as HashMap;
use petgraph::graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph};
use petgraph::Direction;
// use std::collections::HashMap;
use std::iter::FromIterator;

pub fn to_node_list<N: Clone, E>(graph: &DiGraph<N, E>) -> Vec<(usize, N)> {
    graph
        .node_indices()
        .map(|node| {
            let weight = graph.node_weight(node).unwrap().clone();
            (node.index(), weight)
        })
        .collect()
}

pub fn to_edge_list<N, E: Clone>(graph: &DiGraph<N, E>) -> Vec<(usize, usize, usize, E)> {
    graph
        .edge_indices()
        .map(|edge| {
            let weight = graph.edge_weight(edge).unwrap().clone();
            let (source, target) = graph.edge_endpoints(edge).unwrap();
            (edge.index(), source.index(), target.index(), weight)
        })
        .collect()
}

///
/// count the number of (in_degree, out_degree) of nodes.
///
pub fn degree_stats<N, E>(graph: &DiGraph<N, E>) -> HashMap<(usize, usize), usize> {
    let mut h = HashMap::default();
    for node in graph.node_indices() {
        let in_degree = graph.edges_directed(node, Direction::Incoming).count();
        let out_degree = graph.edges_directed(node, Direction::Outgoing).count();
        *h.entry((in_degree, out_degree)).or_insert(0) += 1;
    }
    h
}

/// Remove edges
///
///
pub fn purge_edges_with_mapping<N, E>(
    graph: &mut DiGraph<N, E>,
    edges: &[EdgeIndex],
) -> (
    HashMap<EdgeIndex, Option<EdgeIndex>>,
    HashMap<EdgeIndex, EdgeIndex>,
) {
    let mut from_original: HashMap<EdgeIndex, Option<EdgeIndex>> = HashMap::default();
    let mut to_original = HashMap::default();

    for &edge_remove_original in edges {
        let edge_remove = match from_original.get(&edge_remove_original) {
            None => edge_remove_original,
            Some(&v) => v.expect("remove edge twice"),
        };

        let edge_swap = EdgeIndex::new(graph.edge_count() - 1);
        let edge_swap_original = to_original.get(&edge_swap).copied().unwrap_or(edge_swap);

        // update graph
        graph
            .remove_edge(edge_remove)
            .expect("edge to be removed does not exist");

        if edge_swap == edge_remove {
            from_original.insert(edge_remove_original, None);
            to_original.remove(&edge_remove);
        } else {
            // update from_original and to_original hashmap
            // edge_remove_original is no longer in the graph
            from_original.insert(edge_remove_original, None);
            // swapped edge is now in the removed position
            from_original.insert(edge_swap_original, Some(edge_remove));
            // the last edge in the new graph no longer exists
            to_original.insert(edge_remove, edge_swap_original);
            to_original.remove(&edge_swap);
        }
    }

    (from_original, to_original)
}

fn assert_is_bimap(
    from_original: &HashMap<EdgeIndex, Option<EdgeIndex>>,
    to_original: &HashMap<EdgeIndex, EdgeIndex>,
) {
    for (x, y) in to_original.iter() {
        assert_eq!(from_original.get(y).unwrap().unwrap(), *x);
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ei;

    #[test]
    fn purge_edges_with_mapping_test() {
        let graph: DiGraph<(), usize> = DiGraph::from_edges(&[
            (0, 1, 0),
            (1, 2, 1),
            (2, 3, 2),
            (4, 3, 3),
            (2, 5, 4),
            (3, 7, 5),
            (5, 7, 6),
            (5, 6, 7),
            (6, 5, 8),
        ]);
        assert_eq!(graph.edge_count(), 9);
        assert_eq!(graph.node_count(), 8);

        {
            let mut g = graph.clone();
            let (f, t) = purge_edges_with_mapping(&mut g, &[ei(0)]);
            println!("from_original={:?} to_original={:?}", f, t);
            println!("{:?}", petgraph::dot::Dot::with_config(&g, &[]));
            assert_is_bimap(&f, &t);
            assert_eq!(f, HashMap::from_iter([(ei(0), None), (ei(8), Some(ei(0)))]));
            assert_eq!(t, HashMap::from_iter([(ei(0), ei(8))]));
            assert_eq!(g.edge_count(), 8);
            assert_eq!(g.node_count(), 8);
        }

        {
            let mut g = graph.clone();
            let (f, t) = purge_edges_with_mapping(&mut g, &[ei(3), ei(8), ei(5)]);
            println!("from_original={:?} to_original={:?}", f, t);
            println!("{:?}", petgraph::dot::Dot::with_config(&g, &[]));
            assert_is_bimap(&f, &t);
            assert_eq!(
                f,
                HashMap::from_iter([
                    (ei(3), None),
                    (ei(8), None),
                    (ei(5), None),
                    (ei(7), Some(ei(3))),
                    (ei(6), Some(ei(5))),
                ])
            );
            assert_eq!(t, HashMap::from_iter([(ei(5), ei(6)), (ei(3), ei(7))]));
            assert_eq!(g.edge_count(), 6);
            assert_eq!(g.node_count(), 8);
        }

        {
            let mut g = graph.clone();
            let (f, t) = purge_edges_with_mapping(
                &mut g,
                &[
                    ei(3),
                    ei(8),
                    ei(5),
                    ei(7),
                    ei(0),
                    ei(4),
                    ei(6),
                    ei(1),
                    ei(2),
                ],
            );
            println!("from_original={:?} to_original={:?}", f, t);
            println!("{:?}", petgraph::dot::Dot::with_config(&g, &[]));
            assert_is_bimap(&f, &t);
            assert_eq!(
                f,
                HashMap::from_iter([
                    (ei(0), None),
                    (ei(1), None),
                    (ei(2), None),
                    (ei(3), None),
                    (ei(4), None),
                    (ei(5), None),
                    (ei(6), None),
                    (ei(7), None),
                    (ei(8), None),
                ])
            );
            assert_eq!(t, HashMap::default());
            assert_eq!(g.edge_count(), 0);
            assert_eq!(g.node_count(), 8);
        }
    }
}
