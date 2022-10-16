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
