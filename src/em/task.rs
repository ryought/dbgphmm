//!
//! EM Task and TaskLog definitions
//!
use super::compression::v3::CompressionV3Log;
use super::compression::CompressionLog;
use super::extension::ExtensionLog;
use crate::common::Freq;
use crate::dbg::dbg::{DbgEdge, DbgNode};
use crate::e2e::Dataset;
use itertools::Itertools; // for .join("\n")

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
    /// V3 compression with lambda (1st f64) and zero_penalty (2nd f64)
    ///
    CompressionV3(f64, f64),
    ///
    /// Run extension to specified k-mer length
    ///
    Extension(usize),
}

impl std::fmt::Display for Task {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Task::Compression(depth) => write!(f, "C(d={})", depth),
            Task::CompressionV3(lambda, zero_penalty) => {
                write!(f, "CV3(l={},L0={})", lambda, zero_penalty)
            }
            Task::Extension(k) => write!(f, "E(k={})", k),
        }
    }
}

#[derive(Clone)]
pub enum TaskLog<N: DbgNode, E: DbgEdge> {
    Compression(Vec<CompressionLog<N, E>>),
    CompressionV3(Vec<CompressionV3Log<N, E>>),
    Extension(Vec<ExtensionLog<N, E>>),
}

//
// string conversion
//

impl<N: DbgNode, E: DbgEdge> std::fmt::Display for TaskLog<N, E> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.to_string_with_header(&""))
    }
}

impl<N: DbgNode, E: DbgEdge> TaskLog<N, E> {
    pub fn to_string_with_header(&self, header: &str) -> String {
        match self {
            TaskLog::Compression(logs) => logs
                .iter()
                .enumerate()
                .map(|(step, log)| format!("{}C\t{}\t{}", header, step, log))
                .join("\n"),
            TaskLog::CompressionV3(logs) => logs
                .iter()
                .enumerate()
                .map(|(step, log)| format!("{}CV3\t{}\t{}", header, step, log))
                .join("\n"),
            TaskLog::Extension(logs) => logs
                .iter()
                .enumerate()
                .map(|(step, log)| format!("{}E\t{}\t{}", header, step, log))
                .join("\n"),
        }
    }
    pub fn to_benchmark_string_with_header(&self, dataset: &Dataset, header: &str) -> String {
        match self {
            TaskLog::Compression(logs) => logs
                .iter()
                .enumerate()
                .map(|(step, log)| format!("{}C\t{}\t{}", header, step, log))
                .join("\n"),
            TaskLog::CompressionV3(logs) => logs
                .iter()
                .enumerate()
                .map(|(step, log)| {
                    format!(
                        "{}CV3\t{}\t{}",
                        header,
                        step,
                        log.to_benchmark_string(dataset)
                    )
                })
                .join("\n"),
            TaskLog::Extension(logs) => logs
                .iter()
                .enumerate()
                .map(|(step, log)| {
                    format!(
                        "{}E\t{}\t{}",
                        header,
                        step,
                        log.to_benchmark_string(dataset)
                    )
                })
                .join("\n"),
        }
    }
}
