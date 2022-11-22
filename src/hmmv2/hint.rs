//!
//! Definition of Hint
//!
//!

use petgraph::graph::NodeIndex;

pub struct Hint();

impl Hint {
    ///
    /// Hint for emission sequences `xs` should supports
    ///
    /// `active_nodes(i)` for `xs[i]`
    ///
    pub fn hint_nodes(&self, index: usize) -> &[NodeIndex] {
        unimplemented!();
    }
}
