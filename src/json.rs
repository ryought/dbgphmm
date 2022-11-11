//!
//! json+serde playground
//!
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct User {
    email: usize,
    id: usize,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct List {
    a: Vec<User>,
    #[serde(flatten)]
    b: Vec<User>,
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_flatten() {
        let list = List {
            a: vec![
                User { email: 10, id: 0 },
                User { email: 20, id: 2 },
                User { email: 30, id: 4 },
            ],
            b: vec![
                User { email: 30, id: 4 },
                User { email: 30, id: 4 },
                User { email: 8, id: 5 },
            ],
        };
        let json = serde_json::to_string_pretty(&list).unwrap();
        println!("{}", json);
    }
}
