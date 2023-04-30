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
pub fn euler_circuit_count(graph: &DiGraph<(), usize>) -> f64 {
    let mut graph = graph.clone();
    graph.retain_edges(|g, e| g[e] > 0);
    graph.retain_nodes(|g, v| g.edges_directed(v, Direction::Outgoing).count() > 0);

    println!("cc={}", connected_components(&graph));
    println!("cc={:?}", tarjan_scc(&graph));

    let n = graph.node_count();
    if n == 0 {
        return f64::NEG_INFINITY;
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
    println!("L={}", laplacian);
    let (sign, ln) = laplacian.sln_det().unwrap();
    println!("{} {}", sign, ln);
    // println!("detL={}", sign * ln.exp());

    //
    // PartB:
    //
    let mut ret = if ln == f64::NEG_INFINITY {
        f64::NEG_INFINITY
    } else {
        sign * ln
    };
    for i in graph.node_indices() {
        let c: usize = graph
            .edges_directed(i, Direction::Outgoing)
            .map(|e| e.weight())
            .sum();
        if c > 0 {
            ret += log_factorial(c - 1);

            for e in graph.edges_directed(i, Direction::Outgoing) {
                let c = e.weight();
                ret -= log_factorial(*c);
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

    #[test]
    fn f() {}
}
