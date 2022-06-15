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
}
*/
