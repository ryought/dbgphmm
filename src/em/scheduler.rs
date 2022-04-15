//!
//! EM scheduler
//!

use crate::common::Freq;

///
/// EM two task
///
pub enum Task {
    ///
    /// Compression to the specified depth
    ///
    Compression(Freq),
    ///
    /// Run extension
    ///
    Extension,
}

///
/// Scheduler is a struct that can assign
/// a proper task to the iteration.
///
pub trait Scheduler {
    fn task(&self, iteration: usize) -> Option<Task>;
}

//
// T1
//

///
/// type-1 scheduler
///
pub struct SchedulerType1 {
    /// initial k-mer size
    k_init: usize,
    /// final target k-mer size
    k_target: usize,
    /// true depth
    true_depth: Freq,
}

impl SchedulerType1 {
    pub fn new(k_init: usize, k_target: usize, true_depth: Freq) -> Self {
        SchedulerType1 {
            k_init,
            k_target,
            true_depth,
        }
    }
}

impl Scheduler for SchedulerType1 {
    fn task(&self, iteration: usize) -> Option<Task> {
        unimplemented!();
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn starting() {}
}
