use arrayvec::ArrayVec;

pub fn test() {
    println!("array vec test");
    let mut array = ArrayVec::<u8, 4>::new();
    array.push(1);
    array.push(3);
    array.push(5);
    println!("{:?} {}", array, array.len());
}
