use serde::{Deserialize, Serialize};
use serde_json;

#[derive(Serialize, Deserialize, Debug)]
pub struct DbgStats {
    x: u32,
    y: u32,
}

#[derive(Serialize, Deserialize, Debug, Default)]
pub struct DegreeStats {
    pub in_degs: [u32; 6],
    pub out_degs: [u32; 6],
}

#[derive(Serialize, Deserialize, Debug, Default)]
pub struct CompareResult {
    pub n_shared: u32,
    pub n_self_only: u32,
    pub n_other_only: u32,
}

pub fn test() {
    let point = DbgStats { x: 1, y: 2 };
    let serialized = serde_json::to_string_pretty(&point).unwrap();
    println!("{}", serialized);
}
