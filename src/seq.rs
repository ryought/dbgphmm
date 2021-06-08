#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
enum Base {
    A,
    C,
    G,
    T,
    N,
}

impl Base {
    fn from_u8(n: u8) -> Base {
        match n {
            b'A' => Base::A,
            b'C' => Base::C,
            b'G' => Base::G,
            b'T' => Base::T,
            _ => Base::N,
        }
    }
    fn to_u8(self) -> u8 {
        match self {
            Base::A => b'A',
            Base::C => b'C',
            Base::G => b'G',
            Base::T => b'T',
            Base::N => b'-',
        }
    }
    fn complement(self) -> Base {
        match self {
            Base::A => Base::T,
            Base::C => Base::G,
            Base::G => Base::C,
            Base::T => Base::A,
            Base::N => Base::N,
        }
    }
}

impl std::fmt::Display for Base {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.to_u8() as char)
    }
}

#[derive(Debug)]
struct Seq {
    bases: Vec<Base>,
}

/*
impl Seq {
    fn from_string() -> Seq {}
    fn reverse_complement() -> Seq {
        text.into_iter().rev().map(|a| a.complement()).collect()
    }
}
*/

impl std::fmt::Display for Seq {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for b in self.bases.iter() {
            write!(f, "{}", b.to_u8() as char)?;
        }
        Ok(())
    }
}

trait BaseSeq {
    fn my_length(&self) -> u32;
}

impl BaseSeq for Vec<Base> {
    fn my_length(&self) -> u32 {
        let mut n = 0;
        for b in self.iter() {
            n += 1;
        }
        n
    }
}

pub fn test() {
    println!("seq test");
    let s = Base::A;
    let x = b'T';
    let t = Base::from_u8(x);
    println!("{} {} {}", s, t, t.complement());

    let s1 = Seq {
        bases: vec![Base::A, Base::T, Base::G, Base::T],
    };
    println!("{}", s1);
    println!("my_length = {}", s1.bases.my_length());

    let txt1 = b"ATATATATA"; // &[u8]
    let txt2 = "ATATATTAATT"; // &str
    println!("{} {}", std::str::from_utf8(txt1).unwrap(), txt2);
}

#[derive(Debug)]
struct MySeq(Vec<Base>);
impl std::fmt::Display for MySeq {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for b in self.0.iter() {
            write!(f, "{}", b.to_u8() as char)?;
        }
        Ok(())
    }
}
pub fn test2() {
    let s1 = MySeq(vec![Base::A, Base::T, Base::G, Base::T]);
    println!("seq {}", s1);
}
