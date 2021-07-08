use serde::{Deserialize, Serialize};
use serde_json;

#[derive(Serialize, Deserialize, Debug)]
struct DbgStats {
    x: u32,
    y: u32,
}

pub fn test() {
    let point = DbgStats { x: 1, y: 2 };
    let serialized = serde_json::to_string_pretty(&point).unwrap();
    println!("{}", serialized);
}
