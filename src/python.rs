//!
//! Python bindings
//!
//! Supported features
//!
//! - kmer
//! - styled sequence
//!

use crate::common::collection::StyledSequence;
use crate::kmer::veckmer::VecKmer;
use pyo3::prelude::*;

//
// Dbg
//
use crate::dbg::{dbg::Dbg, impls::SimpleDbg};
#[pyclass]
struct PyDbg(SimpleDbg<VecKmer>);
#[pymethods]
impl PyDbg {
    // #[new]
    // pub fn new(k: usize, seqs: &[&str]) -> Self {
    //     let dbg = Dbg::from_seqs(k, &styled_seqs);
    //     PyDbg(dbg)
    // }
    pub fn k(&self) -> usize {
        self.0.k()
    }
}

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a * a + b * b).to_string())
}

/// A Python module implemented in Rust.
#[pymodule]
fn dbgphmm(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_class::<PyDbg>()?;
    m.add_class::<VecKmer>()?;
    m.add_class::<StyledSequence>()?;
    Ok(())
}
