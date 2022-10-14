// use crate::hoge::{MyTuple, MyVec};
use crate::vector::dense::{DenseFloatStorage, DenseIntStorage};
use pyo3::prelude::*;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a * a + b * b).to_string())
}

/// A Python module implemented in Rust.
#[pymodule]
fn dbgphmm(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_class::<DenseIntStorage>()?;
    m.add_class::<DenseFloatStorage>()?;
    Ok(())
}
