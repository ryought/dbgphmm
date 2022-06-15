//!
//! **Deprecated**
//!
//! Cusom bellman fords
//!
use super::FloatWeight;
use crate::graph::min_mean_weight_cycle::ShortestPaths;
use petgraph::prelude::*;
use petgraph::visit::{VisitMap, Visitable};

//
// tests
//

/*
#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};

    #[test]
    fn shortest_paths_by_edge_01() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            (0, 1, 1.0),
            (1, 2, 4.0),
            (2, 3, 10.0),
            (3, 0, 2.0),
            (3, 4, 2.0),
            (4, 5, 2.0),
            (5, 3, 2.0),
        ]);
        let s = ni(0);
        let p = shortest_paths_by_edge(&g, s, |_, _| true);
        println!("{:?}", p);
        let m = find_minimizer_pair(&g, &p);
        println!("min={:?}", m);
        let (_, e, _) = m.unwrap();
        let path = traceback_preds(&g, e, &p);
        println!("path={:?}", path);
        let cycle = find_cycle(&g, &path);
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some(vec![ni(4), ni(5), ni(3)]));

        // prohibit use of e4
        let p = shortest_paths_by_edge(&g, ni(0), |_, t| t != ei(4));
        println!("{:?}", p);
        let m = find_minimizer_pair(&g, &p);
        println!("min={:?}", m);
        let (_, e, _) = m.unwrap();
        let path = traceback_preds(&g, e, &p);
        println!("path={:?}", path);
        let cycle = find_cycle(&g, &path);
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some(vec![ni(0), ni(1), ni(2), ni(3)]));
    }
    #[test]
    fn mmwc_edge_01() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            (0, 1, 1.0),
            (1, 2, 3.0),
            (2, 0, 1.0),
            (1, 3, 1.0),
            (3, 4, 2.0),
            (4, 5, 1.0),
            (5, 1, 1.0),
        ]);
        let cycle = find_minimum_mean_weight_cycle(&g, ni(0), |_, _| true);
        println!("{:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(3), ni(4), ni(5), ni(1)], 1.25)));
    }
    #[test]
    fn mmwc_edge_02() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            (0, 1, 3.0),
            (0, 1, 1.0),
            (1, 2, 1.0),
            (2, 0, 1.0),
            (1, 3, 1.0),
            (3, 2, 4.0),
        ]);
        let cycle = find_minimum_mean_weight_cycle(&g, ni(0), |_, _| true);
        println!("{:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(0), ni(1), ni(2)], 1.0)));
    }
    #[test]
    fn mmwc_edge_03() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            (0, 1, 1.0),
            (1, 2, 3.0),
            (2, 0, -1.0),
            (1, 3, 2.0),
            (3, 4, 1.0),
            (4, 5, -1.0),
            (5, 6, 2.0),
            (6, 1, 1.0),
            (3, 7, 1.0),
            (7, 4, 2.0),
        ]);

        let cycle = find_minimum_mean_weight_cycle(&g, ni(0), |_, _| true);
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(2), ni(0), ni(1)], 1.0)));

        let cycle = find_minimum_mean_weight_cycle(&g, ni(0), |_, e| e != ei(1));
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(6), ni(1), ni(3), ni(4), ni(5)], 1.0)));
    }
    #[test]
    fn mmwc_edge_04() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            (0, 1, -5.0),
            (1, 2, -5.0),
            (2, 3, -5.0),
            (3, 4, -5.0),
            (4, 0, -5.0),
            (1, 0, -10.0),
            (2, 1, 10.0),
            (3, 2, 10.0),
            (4, 3, 10.0),
            (0, 4, 10.0),
        ]);

        // (1) mmwc among all cycles
        let cycle = find_minimum_mean_weight_cycle(&g, ni(0), |_, _| true);
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(0), ni(1)], -7.5)));

        // (2) restricted mmwc
        let cycle = find_minimum_mean_weight_cycle(&g, ni(0), |e_a, e_b| {
            e_a.index().abs_diff(e_b.index()) != 5
        });
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(2), ni(3), ni(4), ni(0), ni(1)], -5.0)));
    }
}
*/
