pub mod seq {
    pub enum MyDNA {
        A,
        T,
        G,
        C,
        N,
    }
    pub fn complement(x: MyDNA) -> MyDNA {
        match x {
            A => MyDNA::T,
            T => MyDNA::A,
            G => MyDNA::C,
            C => MyDNA::G,
            N => MyDNA::N,
        }
    }
    pub fn show(x: MyDNA) -> char {
        match x {
            A => 'A',
            T => 'T',
            G => 'G',
            C => 'C',
            N => '-',
        }
    }
    pub struct MyDNASeq {}
}

pub mod stat {
    pub struct MyLogProb {}
}

pub mod fuga {
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

pub mod mystring {
    pub fn test1() {
        let s = String::from("hoge");
        println!("{}", s);
    }
}
