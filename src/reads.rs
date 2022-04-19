//!
//! Read/Sequence collection traits
//!
pub fn sequence_to_string(seq: &Bases) -> &str {
    std::str::from_utf8(seq).unwrap()
}

pub type Bases = [u8];

// Empty traits
pub trait Seq: AsRef<Bases> {
    fn to_str(&self) -> &str {
        std::str::from_utf8(self.as_ref()).unwrap()
    }
}
impl<T: AsRef<Bases>> Seq for T {}

fn bar<T>(xs: T)
where
    T: IntoIterator,
    T::Item: Seq,
{
    for x in xs {
        // println!("bar={}", sequence_to_string(x.as_ref()));
        println!("bar={} L={}", x.to_str(), x.as_ref().len());
    }
}

fn foo<'a, I: IntoIterator<Item = &'a u8>>(xs: I) {
    for x in xs {
        println!("{:?}", x);
    }
}

struct Reads(Vec<Vec<u8>>);

impl Reads {
    fn new(reads: Vec<Vec<u8>>) -> Self {
        Reads(reads)
    }
}

impl<'a> IntoIterator for &'a Reads {
    type Item = &'a Vec<u8>;
    type IntoIter = std::slice::Iter<'a, Vec<u8>>;
    fn into_iter(self) -> std::slice::Iter<'a, Vec<u8>> {
        self.0.iter()
    }
}

struct SampledRead {
    bases: Vec<u8>,
    id: usize,
}

impl SampledRead {
    fn new(bases: Vec<u8>, id: usize) -> Self {
        SampledRead { bases, id }
    }
}

struct SampledReads(Vec<SampledRead>);

impl SampledReads {
    fn new(reads: Vec<SampledRead>) -> Self {
        SampledReads(reads)
    }
}

impl<'a> IntoIterator for &'a SampledReads {
    type Item = &'a SampledRead;
    type IntoIter = std::slice::Iter<'a, SampledRead>;
    fn into_iter(self) -> std::slice::Iter<'a, SampledRead> {
        self.0.iter()
    }
}

impl AsRef<Bases> for SampledRead {
    fn as_ref(&self) -> &Bases {
        &self.bases
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn into_iterator_reads() {
        let xs = vec![b"ATCGT", b"TTCGT"];
        bar(&xs);
        bar(&[b"CCCC".to_vec(), b"ATT".to_vec()]);
        bar(&[b"AAAAAAA"]);
        let rs = Reads::new(vec![b"ATCGT".to_vec(), b"TTCGTT".to_vec()]);
        bar(&rs);
        let ss = SampledReads::new(vec![
            SampledRead::new(b"AAAAA".to_vec(), 0),
            SampledRead::new(b"CTCGCTCTC".to_vec(), 1),
        ]);
        bar(&ss);
    }
}
