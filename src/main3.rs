extern crate bio;
use bio::io::fasta;
use bio::stats::LogProb;
use std::io;

// lifetime
#[derive(Debug)]
struct Person<'a> {
    name: &'a str,
    age: u8,
}

// box
struct MyBox {
    b: Box<i32>,
}
impl Drop for MyBox {
    fn drop(&mut self) {
        println!("drop called");
    }
}
fn create_box() {
    let _box: Box<i32> = Box::new(3);
}
fn get_box_1(b: Box<i32>) {
    println!("got box {}", b);
}
fn get_box_2(b: &Box<i32>) {
    println!("borrowed box {}", b);
}

struct Hoge {
    id: String,
    x: u64,
    y: u64,
}

impl std::fmt::Display for Hoge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}:({},{})", self.id, self.x, self.y)
    }
}

#[derive(Debug)]
struct Prob {
    log: f64,
}
use std::ops;
impl ops::Add<Prob> for Prob {
    type Output = Prob;
    fn add(self, _rhs: Prob) -> Prob {
        Prob {
            log: self.log + _rhs.log,
        }
    }
}

// generics
// Copy or Clone
// TODO List all prelude traits
fn largest<T: PartialOrd + Copy>(list: &[T]) -> T {
    let mut largest = list[0];
    for &item in list.iter() {
        if item > largest {
            largest = item;
        }
    }
    largest
}

fn main() {
    /*
    let args: Vec<String> = std::env::args().collect();
    if args.len() == 2 {
        let reader = fasta::Reader::from_file(&std::path::Path::new(&args[1]))
        for record in reader.records() {
            let record = record.unwrap();
            println!("{}, {}", record.id(), record.seq())
        }
    } else {
        println!("please give me a fasta filename");
    }
    */
    let hoge = Hoge {
        id: "hoge22".to_string(),
        x: 0,
        y: 10,
    };
    println!("{}", hoge);

    // let p = LogProb(0.005);
    // let q = LogProb(0.005);
    // println!("{} + {} = {}", p, q, p + q);

    let p = Person {
        name: "hoge",
        age: 10,
    };
    println!("{:?} {}", p, p.age);

    for i in 0u32..1000 {
        // println!("create box {}", i);
        create_box();
    }

    {
        println!("MyBox");
        let x = MyBox { b: Box::new(5) };
        println!("MyBox {}", x.b);
    }

    // mut??
    let b: Box<i32> = Box::new(5);
    get_box_2(&b);
    // get_box_1(b);
    // println!("{}", b);
    let immut_box = Box::new(111u32);
    // *immut_box = 10;
    let mut mut_box = immut_box;
    // let mut_box = 10;
    // mut_box = Box::new(122u32);
    *mut_box = 10;
    println!("{}", mut_box);

    // iterator
    for (i, arg) in std::env::args().enumerate() {
        println!("arg #{} {}", i, arg);
    }
}
