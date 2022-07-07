//!
//! EM algorithms
//!
//! * Compression (deprecated, v1)
//! * CompressionV3
//! * Extension
//!
//! ## Displaying
//!
//! * to_string (std::fmt::Display)
//!     for showing without dataset
//!
//! * to_benchmark_string
//!     for showing with dataset and true genome information
//!
pub mod compression;
pub mod e2e;
pub mod extension;
pub mod scheduler;
pub mod task;
use crate::common::{CopyNum, Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, EdgeCopyNums};
use crate::hmmv2::params::PHMMParams;
use scheduler::Scheduler;
pub use scheduler::SchedulerType1;
pub use task::{Task, TaskLog};

///
/// Do EM inference of Dbg.
/// without iteration end callback
///
pub fn infer<N: DbgNode, E: DbgEdge, S: Scheduler>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    scheduler: &S,
    genome_size: CopyNum,
    max_iter: usize,
) -> (Dbg<N, E>, Vec<TaskLog<N, E>>) {
    infer_with_on_iteration(
        dbg,
        reads,
        params,
        scheduler,
        genome_size,
        max_iter,
        |_, _, _, _| {},
    )
}

///
/// Do EM inference of Dbg.
///
/// with a callback on each end of iteration
/// callback function should take (iteration: usize, Task, &TaskLog, &Dbg).
///
pub fn infer_with_on_iteration<
    N: DbgNode,
    E: DbgEdge,
    S: Scheduler,
    F: Fn(usize, Task, &TaskLog<N, E>, &Dbg<N, E>),
>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    scheduler: &S,
    genome_size: CopyNum,
    max_iter: usize,
    on_iteration: F,
) -> (Dbg<N, E>, Vec<TaskLog<N, E>>) {
    let mut dbg = dbg.clone();
    let mut logs = Vec::new();

    for iteration in 0..scheduler.n_tasks() {
        let task = scheduler.task(iteration).unwrap();
        // run task
        let (dbg_new, log) = match task {
            Task::Compression(depth) => {
                let (dbg_new, log) =
                    compression::v1::compression(&dbg, reads, params, depth, max_iter);
                (dbg_new, TaskLog::Compression(task, log))
            }
            Task::CompressionV3(lambda, zero_penalty) => {
                let (dbg_new, log) = compression::v3::compression(
                    &dbg,
                    reads,
                    params,
                    genome_size,
                    lambda,
                    zero_penalty,
                    max_iter,
                    max_iter,
                );
                (dbg_new, TaskLog::CompressionV3(task, log))
            }
            Task::Extension(k) => {
                let (dbg_new, log) = extension::extension(&dbg, reads, params, max_iter);
                (dbg_new, TaskLog::Extension(task, log))
            }
        };
        // callback
        on_iteration(iteration, task, &log, &dbg_new);

        logs.push(log);
        dbg = dbg_new;
    }

    (dbg, logs)
}
