//!
//! Greedy search of posterior probability
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase, EdgeCopyNums, NodeCopyNums};
use crate::e2e::Dataset;
use crate::greedy::{GreedyInstance, GreedySearcher};
use crate::hmmv2::params::PHMMParams;

#[derive(Clone, Hash, PartialEq, Eq)]
struct DbgCopyNumsInstance(NodeCopyNums);
impl DbgCopyNumsInstance {
    fn new(cn: NodeCopyNums) -> Self {
        DbgCopyNumsInstance(cn)
    }
    fn copy_nums(&self) -> &NodeCopyNums {
        &self.0
    }
}
impl GreedyInstance for DbgCopyNumsInstance {}

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    pub fn search_posterior_once(&self, dataset: &Dataset) {
        let instance_init = DbgCopyNumsInstance::new(self.to_node_copy_nums());
        let mut searcher = GreedySearcher::new(
            instance_init,
            |instance| {
                let mut dbg = self.clone();
                dbg.set_node_copy_nums(instance.copy_nums());
                dbg.to_full_prob(dataset.params(), dataset.reads())
            },
            |instance| {
                let mut dbg = self.clone();
                dbg.set_node_copy_nums(instance.copy_nums());
                dbg.neighbor_copy_nums_fast_compact(100)
                    .into_iter()
                    .map(|copy_nums| DbgCopyNumsInstance::new(copy_nums))
                    .collect()
            },
        );
        let n_new_neighbors = searcher.search_once();
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::e2e::generate_small_tandem_repeat;

    #[test]
    fn aaa() {}
}
