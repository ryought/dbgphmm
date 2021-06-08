mod hoge {
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
    let y = crate::hoge::foo::hoge(20);
    println!("y {}", y);
}
