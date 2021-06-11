pub mod hmm;
pub mod hoge; // load file src/hoge.rs
pub mod kmer;
pub mod linkedlist;
pub mod my_vec;
pub mod prob;
pub mod sleeper;
pub mod test_struct;
pub mod vec_of_vec;

#[macro_use]
extern crate approx;

mod fuga {
    pub mod foo {
        pub fn hoge(x: i32) -> i32 {
            x + 100
        }
    }
    pub fn largest<T: PartialOrd + Copy>(list: &[T]) -> T {
        let mut largest = list[0];
        for &item in list.iter() {
            if item > largest {
                largest = item;
            }
        }
        largest
    }
}

pub fn test() {
    let y = crate::fuga::foo::hoge(20);
    println!("y {}", y);
}
