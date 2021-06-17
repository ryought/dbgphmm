#[derive(Debug)]
struct Hoge {
    hoge: Vec<Vec<u8>>,
}

impl Hoge {
    fn new() -> Hoge {
        let mut v = Vec::new();
        for _ in 0..10 {
            let w = vec![1, 2, 3, 4, 5];
            v.push(w);
        }
        Hoge { hoge: v }
    }
    // if it returns the reference to the member, lifetime is needed
    fn get_ref<'a>(&'a self, index: usize) -> &'a [u8] {
        &self.hoge[index]
    }
    // if it returns the iterator to the reference of the member
    fn get_iter(&self, index: usize) -> impl Iterator<Item = u8> + '_ {
        self.hoge[index].iter().map(|x| x + 1)
    }
    /*
    fn get_ref_iter(&self) -> impl Iterator<Item = &Vec<u8>> + '_ {
        self.hoge.iter()
    }
    */
}

/*
trait Fuga<'a> {
    type ItemIterator: Iterator<Item = &'a u8>;
    fn get_ref_iter(&'a self) -> impl Self::ItemIterator;
}

impl<'a> Fuga<'a> for Hoge {
    type ItemIterator;
    fn get_ref_iter(&'a self) -> impl Self::ItemIterator {
        self.hoge.iter()
    }
}

pub fn test() {
    let h = Hoge::new();
    println!("{:?}", h);
    println!("{:?}", h.get_ref(1));
    // h.get_ref(1);
    for i in h.get_iter(1) {
        println!("it: {}", i);
    }
    for i in h.get_ref_iter() {
        println!("it: {:?}", i);
        // i[0] = 10;
    }
}
*/
