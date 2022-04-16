//!
//! EM scheduler
//!

use crate::common::Freq;

///
/// EM two task
///
#[derive(Clone, Debug, PartialEq)]
pub enum Task {
    ///
    /// Compression to the specified depth
    ///
    Compression(Freq),
    ///
    /// Run extension to specified k-mer length
    ///
    Extension(usize),
}

///
/// Scheduler is a struct that can assign
/// a proper task to the iteration.
///
pub trait Scheduler {
    fn n_tasks(&self) -> usize;
    fn task(&self, iteration: usize) -> Option<Task>;
}

//
// T1
//

///
/// type-1 scheduler
///
#[derive(Clone, Debug)]
pub struct SchedulerType1 {
    /// initial k-mer size
    k_init: usize,
    /// final target k-mer size
    k_target: usize,
    /// true depth
    true_depth: Freq,
    //
    // internal variables
    //
    n_extension: usize,
    n_compression: usize,
    n_iteration: usize,
    compression_interval: usize,
}

impl SchedulerType1 {
    pub fn new(k_init: usize, k_target: usize, true_depth: Freq) -> Self {
        assert!(k_target >= k_init);
        assert!(true_depth > 1.0);

        // (1) compute n_extension and n_compression.
        // how many times extension and compression should be executed?
        let n_extension = k_target - k_init;
        let n_compression = (true_depth - 1.0).ceil() as usize;
        let n_iteration = n_extension + n_compression;
        let compression_interval = if n_compression > 0 {
            n_extension / n_compression
        } else {
            0
        } + 1;

        SchedulerType1 {
            k_init,
            k_target,
            true_depth,
            n_extension,
            n_compression,
            n_iteration,
            compression_interval,
        }
    }
}

impl Scheduler for SchedulerType1 {
    fn n_tasks(&self) -> usize {
        self.n_iteration
    }
    fn task(&self, iteration: usize) -> Option<Task> {
        let stage = iteration / self.compression_interval;
        let step_in_stage = iteration % self.compression_interval;
        let is_final_stage =
            (self.n_iteration - stage * self.compression_interval) < self.compression_interval;

        if iteration < self.n_iteration {
            if !is_final_stage {
                if step_in_stage == 0 {
                    // do compression
                    let depth = (2.0 + stage as Freq).min(self.true_depth);
                    Some(Task::Compression(depth))
                } else {
                    // do extension
                    let k = self.k_init + (iteration - stage);
                    Some(Task::Extension(k))
                }
            } else {
                let k = self.k_init + (iteration - stage) + 1;
                Some(Task::Extension(k))
            }
        } else {
            None
        }
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn scheduler_type1() {
        let s = SchedulerType1::new(16, 19, 2.2);
        let tasks: Vec<Task> = (0..s.n_tasks()).map(|i| s.task(i).unwrap()).collect();
        println!("{:?}", tasks);
        assert_eq!(
            tasks,
            vec![
                Task::Compression(2.0),
                Task::Extension(17),
                Task::Compression(2.2),
                Task::Extension(18),
                Task::Extension(19)
            ]
        );

        let s = SchedulerType1::new(16, 32, 5.0);
        let tasks: Vec<Task> = (0..s.n_tasks()).map(|i| s.task(i).unwrap()).collect();
        println!("{:?}", tasks);
        assert_eq!(
            tasks,
            vec![
                Task::Compression(2.0),
                Task::Extension(17),
                Task::Extension(18),
                Task::Extension(19),
                Task::Extension(20),
                Task::Compression(3.0),
                Task::Extension(21),
                Task::Extension(22),
                Task::Extension(23),
                Task::Extension(24),
                Task::Compression(4.0),
                Task::Extension(25),
                Task::Extension(26),
                Task::Extension(27),
                Task::Extension(28),
                Task::Compression(5.0),
                Task::Extension(29),
                Task::Extension(30),
                Task::Extension(31),
                Task::Extension(32)
            ]
        );

        let s = SchedulerType1::new(16, 16, 5.4);
        let tasks: Vec<Task> = (0..s.n_tasks()).map(|i| s.task(i).unwrap()).collect();
        println!("{:?}", tasks);
        assert_eq!(
            tasks,
            vec![
                Task::Compression(2.0),
                Task::Compression(3.0),
                Task::Compression(4.0),
                Task::Compression(5.0),
                Task::Compression(5.4)
            ]
        );
    }
}
