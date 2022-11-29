//!
//! Hint information
//!
use crate::graph::active_nodes::ActiveNodes;

#[derive(Clone, Debug, PartialEq)]
pub struct Hint(Vec<ActiveNodes>);

impl Hint {
    pub fn new(active_nodes_vec: Vec<ActiveNodes>) -> Self {
        Hint(active_nodes_vec)
    }
    pub fn active_nodes(&self, index: usize) -> &ActiveNodes {
        &self.0[index]
    }
    pub fn len(&self) -> usize {
        self.0.len()
    }
}
