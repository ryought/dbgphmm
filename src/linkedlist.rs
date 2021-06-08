/// size of enum
/// sizeof(tag) + max(sizeof(attributes))
#[derive(Debug)]
enum List {
    Cons(i32, Box<List>),
    Nil,
}

use List::{Cons, Nil};

pub fn test1() {
    let l = Cons(1, Box::new(Cons(2, Box::new(Cons(3, Box::new(Nil))))));
    println!("{:?}", l);
}

struct MyBox<T>(T);

impl<T> MyBox<T> {
    fn new(x: T) -> MyBox<T> {
        MyBox(x)
    }
}

pub fn test2() {
    let x = 5;
    let y = &x;
    let y = MyBox::new(x);
    let y = Box::new(x);
    println!("{:?} {:?} {}", x, y, x == *y);
}
