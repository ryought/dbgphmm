// iterator implementation test
struct Counter {
    count: u32,
}
impl Counter {
    fn new() -> Counter {
        Counter { count: 0 }
    }
}
impl Iterator for Counter {
    type Item = u32;
    fn next(&mut self) -> Option<Self::Item> {
        self.count += 1;
        if self.count < 6 {
            Some(self.count)
        } else {
            None
        }
    }
}

fn test2(vec: &[u8]) -> Vec<u8> {
    vec.iter().map(|x| x + 1).collect()
}

fn test() {
    let mut counter = Counter::new();
    for x in counter {
        println!("{:?}", x);
    }

    let vec: Vec<u8> = vec![5, 8, 9];
    let ids: Vec<usize> = (0..vec.len()).collect();
    println!("{:?}, {:?}", vec, ids);
}
