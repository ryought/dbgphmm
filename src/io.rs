pub mod cytoscape;
pub mod fasta;
pub mod json;

use std::fs::File;
use std::io::prelude::*;

///
/// write string into a file
///
pub fn write_string(filename: &str, string: &str) -> std::io::Result<()> {
    let mut file = File::create(&filename)?;
    file.write_all(string.as_bytes())?;
    Ok(())
}
