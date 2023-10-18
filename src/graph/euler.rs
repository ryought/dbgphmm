//!
//! Calculating the number of euler circuits
//!
use crate::utils::log_factorial;
use ndarray::prelude::*;
use ndarray_linalg::solve::Determinant;
use petgraph::algo::connected_components;
use petgraph::algo::tarjan_scc;
use petgraph::graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph};
use petgraph::Direction;

///
///
///
pub fn subgraph<N: Clone, E: Clone>(graph: &DiGraph<N, E>, nodes: &[NodeIndex]) -> DiGraph<N, E> {
    let mut graph = graph.clone();
    graph.retain_nodes(|_, node| nodes.contains(&node));
    graph
}

pub fn euler_circuit_count_in_connected(graph: &DiGraph<(), usize>) -> f64 {
    // println!("[euler] {:?}", petgraph::dot::Dot::with_config(&graph, &[]));
    let n = graph.node_count();
    if n == 0 {
        return 0.0;
    }

    //
    // PartA: create laplacian matrix L
    //
    let mut laplacian: Array2<f64> = Array::zeros((n, n));
    // (1) diag (degree) matrix
    // L[i,i] += (total copy numbers of node i)
    for i in graph.node_indices() {
        let c: usize = graph
            .edges_directed(i, Direction::Outgoing)
            .map(|e| e.weight())
            .sum();
        // println!("i={} c={}", i.index(), c);
        laplacian[[i.index(), i.index()]] = c as f64;
    }
    // (2) subtract adjacency matrix a_ij
    // L[i,j] -= (total copy numbers of edges i->j)
    for i in graph.node_indices() {
        for j in graph.node_indices() {
            let c: usize = graph.edges_connecting(i, j).map(|e| e.weight()).sum();
            // println!("i={} j={} c={}", i.index(), j.index(), c);
            laplacian[[i.index(), j.index()]] -= c as f64;
        }
    }
    // add +1
    // starting point is arbitrary
    laplacian[[0, 0]] += 1.0;
    // println!("L={}", laplacian);
    let (sign, ln) = laplacian.sln_det().unwrap();
    // println!("{} {}", sign, ln);
    // println!("detL={}", sign * ln.exp());

    //
    // PartB:
    //
    let mut count = if ln == f64::NEG_INFINITY {
        0.0
    } else {
        sign * ln
    };
    for i in graph.node_indices() {
        let c: usize = graph
            .edges_directed(i, Direction::Outgoing)
            .map(|e| e.weight())
            .sum();
        if c > 0 {
            count += log_factorial(c - 1);

            for e in graph.edges_directed(i, Direction::Outgoing) {
                let c = e.weight();
                count -= log_factorial(*c);
            }
        }
    }

    count
}

///
/// Count the number of Eulerian circuits in the graph (whose edge has its multiplicity) (in log
/// space)
///
/// if `allow_multiple_component=true`, the final count will be the product of count of each components.
/// otherwise, the resulting count is log(0) if the graph is not strongly connected. (because there
/// is no Eulerian circuit)
///
pub fn euler_circuit_count(graph: &DiGraph<(), usize>, allow_multiple_component: bool) -> f64 {
    let mut graph = graph.clone();
    // remove zero edges
    graph.retain_edges(|g, e| g[e] > 0);
    // remove isolated nodes (= no outgoing edge)
    graph.retain_nodes(|g, v| g.edges_directed(v, Direction::Outgoing).count() > 0);

    // println!("cc={}", connected_components(&graph));
    // println!("cc={:?}", tarjan_scc(&graph));

    let mut ret = 0.0;
    let n = graph.node_count();
    if n == 0 {
        return f64::NEG_INFINITY;
    }

    // for each strongly connected component..
    if allow_multiple_component {
        for component in tarjan_scc(&graph) {
            let h = subgraph(&graph, &component);
            let count = euler_circuit_count_in_connected(&h);
            // println!("compont={:?} count={}", component, count.exp());
            ret += count;
        }
    } else {
        if tarjan_scc(&graph).len() > 1 {
            // if there are separate component, it has no euler circuit
            return f64::NEG_INFINITY;
        } else {
            ret = euler_circuit_count_in_connected(&graph);
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

    #[test]
    fn subgraph_test() {
        let mut g = DiGraph::new();
        let a = g.add_node(());
        let b = g.add_node(());
        let c = g.add_node(());
        g.add_edge(a, b, 1);
        g.add_edge(a, c, 2);
        g.add_edge(c, b, 3);
        g.add_edge(b, c, 4);
        g.add_edge(b, c, 5);
        let h = subgraph(&g, &[b, c]);
        println!("{:?}", petgraph::dot::Dot::with_config(&h, &[]));
    }

    #[test]
    fn n_euler() {
        let assert_euler_approx_eq = |graph: &DiGraph<(), usize>, n0: usize, n1: usize| {
            assert!((euler_circuit_count(&graph, false).exp() - n0 as f64).abs() < 0.001);
            assert!((euler_circuit_count(&graph, true).exp() - n1 as f64).abs() < 0.001);
        };

        // self loop
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[(0, 0, 1)]);
        assert_euler_approx_eq(&g, 1, 1);

        // single edge has no euler circuit
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[(0, 1, 2)]);
        assert_euler_approx_eq(&g, 0, 0);

        // loop
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[(0, 1, 1), (1, 0, 1)]);
        assert_euler_approx_eq(&g, 1, 1);

        // three euler circuit
        // - A,X,X,B
        // - A,X,B,X
        // - A,B,X,X
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[
            (0, 0, 1), // A
            (0, 0, 1), // B
            (0, 0, 2), // X
        ]);
        assert_euler_approx_eq(&g, 3, 3);

        // two bubbles with two euler circuit
        // - A1 B C1 D A2 B C2 D
        // - A1 B C2 D A2 B C1 D
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[
            (0, 1, 1), // A1
            (0, 1, 1), // A2
            (1, 2, 2), // B
            (2, 3, 1), // C1
            (2, 3, 1), // C1
            (3, 0, 2), // D
        ]);
        assert_euler_approx_eq(&g, 2, 2);

        // separate components
        // if allowing separate components, it has euler circuit A and B
        // otherwise it has no circuit
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[
            (0, 0, 1), // A
            (1, 1, 1), // B
        ]);
        assert_euler_approx_eq(&g, 0, 1);

        println!(
            "{} {}",
            euler_circuit_count(&g, false).exp(),
            euler_circuit_count(&g, true).exp(),
        );
    }
}
