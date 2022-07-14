//!
//! EM scheduler
//!

use super::task::Task;
use crate::common::Freq;

///
/// Scheduler is a struct that can assign
/// a proper task to the iteration.
///
pub trait Scheduler {
    fn n_tasks(&self) -> usize;
    fn task(&self, iteration: usize) -> Option<Task>;
}

//
// common utility functions
//

///
/// n items [x0, ...., x1]
/// 0-th item is x0
/// n-1-th item is x1
///
/// then return i-th item
///
fn interpolate_linear(x0: f64, x1: f64, n: usize, i: usize) -> f64 {
    assert!(n >= 1);
    assert!(i < n);
    if n == 1 {
        x0
    } else {
        x0 + (x1 - x0) * i as f64 / (n - 1) as f64
    }
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
    //
    // for v1 compression
    //
    /// true depth
    true_depth: Freq,
    //
    // for v3 compression
    //
    use_v3: bool,
    lambda_init: f64,
    lambda_target: f64,
    zero_penalty_init: f64,
    zero_penalty_target: f64,
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
            // v3 parameters
            use_v3: false,
            lambda_init: 0.0,
            lambda_target: 0.0,
            zero_penalty_init: 0.0,
            zero_penalty_target: 0.0,
        }
    }
    pub fn new_v3(
        k_init: usize,
        k_target: usize,
        n_compression: usize,
        lambda_init: f64,
        lambda_target: f64,
        zero_penalty_init: f64,
        zero_penalty_target: f64,
    ) -> Self {
        assert!(k_target >= k_init);

        // (1) compute n_extension and n_compression.
        // how many times extension and compression should be executed?
        let n_extension = k_target - k_init;
        let n_iteration = n_extension + n_compression;
        let compression_interval = if n_compression > 0 {
            n_extension / n_compression
        } else {
            0
        } + 1;

        SchedulerType1 {
            k_init,
            k_target,
            n_extension,
            n_compression,
            n_iteration,
            compression_interval,
            // v3 parameter
            use_v3: true,
            lambda_init,
            lambda_target,
            zero_penalty_init,
            zero_penalty_target,
            // v1 parameter
            true_depth: 1.0, // temp
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
        let is_final_stage = iteration >= self.n_compression * self.compression_interval;

        if iteration < self.n_iteration {
            if !is_final_stage {
                if step_in_stage == 0 {
                    // do compression
                    if self.use_v3 {
                        let lambda = interpolate_linear(
                            self.lambda_init,
                            self.lambda_target,
                            self.n_compression,
                            stage,
                        );
                        let zero_penalty = interpolate_linear(
                            self.zero_penalty_init,
                            self.zero_penalty_target,
                            self.n_compression,
                            stage,
                        );
                        Some(Task::CompressionV3(lambda, zero_penalty))
                    } else {
                        let depth = (2.0 + stage as Freq).min(self.true_depth);
                        Some(Task::Compression(depth))
                    }
                } else {
                    // do extension
                    let k = self.k_init + (iteration - stage);
                    Some(Task::Extension(k))
                }
            } else {
                let k = self.k_target - (self.n_iteration - iteration - 1);
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
    fn scheduler_interpolate() {
        println!("{}", interpolate_linear(0.0, 10.0, 11, 5));
        assert_eq!(interpolate_linear(0.0, 10.0, 11, 5), 5.0);
        assert_eq!(interpolate_linear(0.0, 10.0, 11, 0), 0.0);
        assert_eq!(interpolate_linear(0.0, 10.0, 11, 10), 10.0);
        assert_eq!(interpolate_linear(0.0, 10.0, 1, 0), 0.0);
    }
    #[test]
    fn scheduler_type1_v1() {
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

        let s = SchedulerType1::new(8, 50, 10.0);
        let tasks: Vec<Task> = (0..s.n_tasks()).map(|i| s.task(i).unwrap()).collect();
        println!("{:?}", tasks);
        assert_eq!(
            tasks,
            vec![
                Task::Compression(2.0),
                Task::Extension(9),
                Task::Extension(10),
                Task::Extension(11),
                Task::Extension(12),
                Task::Compression(3.0),
                Task::Extension(13),
                Task::Extension(14),
                Task::Extension(15),
                Task::Extension(16),
                Task::Compression(4.0),
                Task::Extension(17),
                Task::Extension(18),
                Task::Extension(19),
                Task::Extension(20),
                Task::Compression(5.0),
                Task::Extension(21),
                Task::Extension(22),
                Task::Extension(23),
                Task::Extension(24),
                Task::Compression(6.0),
                Task::Extension(25),
                Task::Extension(26),
                Task::Extension(27),
                Task::Extension(28),
                Task::Compression(7.0),
                Task::Extension(29),
                Task::Extension(30),
                Task::Extension(31),
                Task::Extension(32),
                Task::Compression(8.0),
                Task::Extension(33),
                Task::Extension(34),
                Task::Extension(35),
                Task::Extension(36),
                Task::Compression(9.0),
                Task::Extension(37),
                Task::Extension(38),
                Task::Extension(39),
                Task::Extension(40),
                Task::Compression(10.0),
                Task::Extension(41),
                Task::Extension(42),
                Task::Extension(43),
                Task::Extension(44),
                Task::Extension(45),
                Task::Extension(46),
                Task::Extension(47),
                Task::Extension(48),
                Task::Extension(49),
                Task::Extension(50)
            ]
        );
    }
    #[test]
    fn scheduler_type1_v3() {
        let s = SchedulerType1::new_v3(16, 16, 10, -10., -100., -1., -10.);
        let tasks: Vec<Task> = (0..s.n_tasks()).map(|i| s.task(i).unwrap()).collect();
        println!("{:?}", tasks);
        assert_eq!(
            tasks,
            vec![
                Task::CompressionV3(-10.0, -1.0),
                Task::CompressionV3(-20.0, -2.0),
                Task::CompressionV3(-30.0, -3.0),
                Task::CompressionV3(-40.0, -4.0),
                Task::CompressionV3(-50.0, -5.0),
                Task::CompressionV3(-60.0, -6.0),
                Task::CompressionV3(-70.0, -7.0),
                Task::CompressionV3(-80.0, -8.0),
                Task::CompressionV3(-90.0, -9.0),
                Task::CompressionV3(-100.0, -10.0),
            ]
        );

        let s = SchedulerType1::new_v3(16, 20, 1, -10., -10., -1., -1.);
        let tasks: Vec<Task> = (0..s.n_tasks()).map(|i| s.task(i).unwrap()).collect();
        println!("{:?}", tasks);
        assert_eq!(
            tasks,
            vec![
                Task::CompressionV3(-10.0, -1.0),
                Task::Extension(17),
                Task::Extension(18),
                Task::Extension(19),
                Task::Extension(20),
            ]
        );
    }
}
