//!
//! Python bindings
//!
//! Supported features
//!
//! - kmer
//! - styled sequence
//!

use crate::common::collection::StyledSequence;
use crate::common::CopyNum;
use crate::hmmv2::params::PHMMParams;
use crate::kmer::veckmer::VecKmer;
use pyo3::prelude::*;

//
// Dbg
//
use crate::dbg::{
    dbg::{Dbg, DbgEdge, DbgNode, DbgNodeBase},
    impls::SimpleDbg,
};
#[pyclass]
struct PyDbg(SimpleDbg<VecKmer>);
impl PyDbg {
    fn dbg(&self) -> &SimpleDbg<VecKmer> {
        &self.0
    }
}
#[pymethods]
impl PyDbg {
    #[new]
    fn __new__(k: usize, seqs: Vec<StyledSequence>) -> Self {
        let dbg = Dbg::from_seqs(k, &seqs);
        PyDbg(dbg)
    }
    fn __repr__(&self) -> String {
        self.dbg().to_string()
    }
    //
    // getters
    //
    #[getter]
    fn k(&self) -> usize {
        self.dbg().k()
    }
    #[getter]
    fn n_nodes(&self) -> usize {
        self.dbg().n_nodes()
    }
    #[getter]
    fn n_edges(&self) -> usize {
        self.dbg().n_edges()
    }
    #[getter]
    fn genome_size(&self) -> usize {
        self.dbg().genome_size()
    }
    //
    // graph
    //
    /// get a list of all nodes as (kmer string, copy number)
    fn nodes(&self) -> Vec<(String, CopyNum)> {
        self.dbg()
            .nodes()
            .map(|(_, w)| (w.kmer().to_string(), w.copy_num()))
            .collect()
    }
    /// get a list of all edges as (source node index, target node index, copy number)
    fn edges(&self) -> Vec<(usize, usize, Option<CopyNum>)> {
        self.dbg()
            .edges()
            .map(|(_, s, t, w)| (s.index(), t.index(), w.copy_num()))
            .collect()
    }
    //
    // PHMM related
    //
    /// get NodeFreqs of reads as PHMM
    /// returned value: `(vec of expected usage of each node, log full prob)`
    fn to_node_freqs(&self, param: PHMMParams, reads: Vec<StyledSequence>) -> (Vec<f64>, f64) {
        let phmm = self.dbg().to_phmm(param);
        let (node_freqs, p) = phmm.to_node_freqs_parallel(&reads);
        (node_freqs.to_inner_vec(), p.to_log_value())
    }
    /// get EdgeFreqs and InitFreqs of reads as PHMM
    /// returned value:
    /// ```text
    /// (
    ///     vec of expected usage of each edge/transition,
    ///     vec of expected usage of each node as a first step,
    ///     log full prob
    /// )
    /// ```
    fn to_edge_and_init_freqs(
        &self,
        param: PHMMParams,
        reads: Vec<StyledSequence>,
    ) -> (Vec<f64>, Vec<f64>, f64) {
        let phmm = self.dbg().to_phmm(param);
        let (edge_freqs, init_freqs, p) = phmm.to_edge_and_init_freqs_parallel(&reads);
        (
            edge_freqs.to_inner_vec(),
            init_freqs.to_inner_vec(),
            p.to_log_value(),
        )
    }
    ///
    /// `prob_matrix[i, {0,1,2}, j]`
    /// = (`i`-th base is emitted from `{0:Match, 1:Ins, 2:Del}` state of `j`-th node)
    ///
    fn to_prob_matrix(
        &self,
        param: PHMMParams,
        read: StyledSequence,
    ) -> Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> {
        let phmm = self.dbg().to_phmm(param);
        let o = phmm.run(read.as_ref());
        o.to_naive_emit_probs()
    }
}

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a * a + b * b).to_string())
}

/// string concat test
#[pyfunction]
fn string_concat(strings: Vec<String>) -> String {
    let mut ret = String::new();
    for string in strings {
        ret += &string;
    }
    ret
}

/// A Python module implemented in Rust.
#[pymodule]
fn dbgphmm(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(string_concat, m)?)?;
    m.add_class::<PyDbg>()?;
    m.add_class::<VecKmer>()?;
    m.add_class::<StyledSequence>()?;
    m.add_class::<PHMMParams>()?;
    Ok(())
}
