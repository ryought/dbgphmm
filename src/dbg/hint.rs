//!
//! Generate Hint information for efficient PHMM calculation on Dbg
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::common::{ei, ni, CopyNum, Freq, ReadCollection, Seq};
use crate::dbg::draft::MAX_COPY_NUM_OF_EDGE;
use crate::dbg::edge_centric::compact::compacted_flow_into_original_flow;
use crate::hmmv2::hint::Hint;
use crate::hmmv2::params::PHMMParams;
use crate::min_flow::{min_cost_flow_convex_fast, total_cost, ConstCost, Cost, FlowEdge};
use crate::vector::graph::flow_to_edgevec;

///
/// Uniform edge type for `Dbg::uniform_copy_nums`
///
/// MinCostFlow edge
/// * [demand,capacity] = [1,MAX]
/// * cost_per_unit = 1
///
#[derive(Clone, Debug)]
struct Uniform;
impl FlowEdge<usize> for Uniform {
    fn demand(&self) -> usize {
        1
    }
    fn capacity(&self) -> usize {
        MAX_COPY_NUM_OF_EDGE
    }
}
impl ConstCost for Uniform {
    fn cost(&self) -> f64 {
        1.0
    }
}

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    /// Generate hints
    ///
    /// ## TODO
    /// * Cloning reads is redundant. It should return Vec<(&S, Hint)>.
    ///
    pub fn generate_hints<'a, S: Seq>(
        &self,
        reads: ReadCollection<S>,
        params: PHMMParams,
    ) -> ReadCollection<S> {
        let mut dbg = self.clone();
        let copy_nums = dbg.uniform_copy_nums();
        dbg.set_node_copy_nums(&copy_nums);
        let phmm = dbg.to_phmm(params);

        phmm.append_hints(reads, false, false)
    }
    ///
    /// Find copy nums in which all edges have non-zero copynum
    /// To be used when creating hints for reads.
    ///
    /// Solves MinCostFlow on compacted edge-centric dbg with each edge have `Uniform` edge
    /// ([l,u]=[1,MAX], c=1).
    ///
    /// ## TODO
    /// * all-one copy_nums is also ok? Flow-constraint will be broken, but it seems not to cause
    /// any problems.
    ///
    fn uniform_copy_nums(&self) -> NodeCopyNums {
        let graph = self.to_compact_edbg_graph();
        let flow_network = graph.map(|_, _| (), |_, _| Uniform);
        min_cost_flow_convex_fast(&flow_network)
            .map(|flow| {
                let flow_in_original =
                    compacted_flow_into_original_flow(self.n_nodes(), &graph, &flow);
                flow_to_edgevec(flow_in_original).switch_index()
            })
            .unwrap()
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks;
    use crate::e2e;
    use crate::utils::timer;

    #[test]
    fn uniform_copy_nums_intersection_small() {
        let dbg = mocks::mock_intersection_small();
        let uniform_copy_nums = dbg.uniform_copy_nums();
        println!("{}", uniform_copy_nums);
        assert_eq!(uniform_copy_nums, NodeCopyNums::new(dbg.n_nodes(), 1));
    }

    #[test]
    fn dbg_hint_simple_genome() {
        let exp = e2e::generate_simple_genome_mock();
        let dbg = exp.dbg_draft.clone().unwrap();
        println!("n_reads={}", exp.reads().len());
        let mut param = exp.phmm_params;
        param.n_active_nodes = 10;

        let reads_with_hints = dbg.generate_hints(exp.reads().clone(), param);

        // run for each reads manually
        let phmm = dbg.to_phmm(param);
        for (read, hint) in reads_with_hints.iter_with_hint() {
            let r1 = phmm.forward(read.as_ref());
            let r2 = phmm.forward_with_hint(read.as_ref(), hint);
            let p1 = r1.full_prob();
            let p2 = r2.full_prob();
            println!("p(dense)={} p(hint)={} diff={}", p1, p2, p1.log_diff(p2));
            assert!(p1.log_diff(p2) < 1.0);
        }

        // run for all reads
        let (p1, t1) = timer(|| dbg.to_full_prob(param, exp.reads()));
        let (p2, t2) = timer(|| dbg.to_full_prob_with_hint(param, &reads_with_hints));
        println!(
            "p_tot(dense)={} p_tot(hint)={} diff={}",
            p1,
            p2,
            p1.log_diff(p2)
        );
        println!("t1={} t2={}", t1, t2);
        assert!(p1.log_diff(p2) < 0.1);
    }
}
