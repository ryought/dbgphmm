pub mod compression;
pub mod e2e;
pub mod extension;
pub mod scheduler;
use crate::common::{Freq, Reads};
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, EdgeCopyNums};
use crate::hmmv2::params::PHMMParams;
pub use scheduler::SchedulerType1;
use scheduler::{Scheduler, Task};

///
/// Do EM inference of Dbg.
///
pub fn infer<N: DbgNode, E: DbgEdge, S: Scheduler>(
    dbg: &Dbg<N, E>,
    reads: &Reads,
    params: &PHMMParams,
    scheduler: &S,
    max_iter: usize,
) -> Dbg<N, E> {
    let mut dbg = dbg.clone();

    for iteration in 0..scheduler.n_tasks() {
        let task = scheduler.task(iteration).unwrap();
        match task {
            Task::Compression(depth) => {
                (dbg, _) = compression::compression(&dbg, reads, params, depth, max_iter);
            }
            Task::Extension(k) => {
                (dbg, _) = extension::extension(&dbg, reads, params, max_iter);
            }
        }
    }

    dbg
}
