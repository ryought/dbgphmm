pub mod compression;
pub mod e2e;
pub mod extension;
pub mod scheduler;
use crate::common::{CopyNum, Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, EdgeCopyNums};
use crate::hmmv2::params::PHMMParams;
pub use scheduler::SchedulerType1;
use scheduler::{Scheduler, Task, TaskLog};

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
        |_, _, _| {},
    )
}

///
/// Do EM inference of Dbg.
///
pub fn infer_with_on_iteration<
    N: DbgNode,
    E: DbgEdge,
    S: Scheduler,
    F: Fn(Task, &[TaskLog<N, E>], &Dbg<N, E>),
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
        match task {
            Task::Compression(depth) => {
                let (dbg_new, log) = compression::compression(&dbg, reads, params, depth, max_iter);
                logs.push(TaskLog::Compression(log));
                dbg = dbg_new;
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
                logs.push(TaskLog::CompressionV3(log));
                dbg = dbg_new;
            }
            Task::Extension(k) => {
                let (dbg_new, log) = extension::extension(&dbg, reads, params, max_iter);
                logs.push(TaskLog::Extension(log));
                dbg = dbg_new;
            }
        }

        on_iteration(task, &logs, &dbg);
    }

    (dbg, logs)
}
