use bio::io::fasta;
use dbgphmm::common::sequence_to_string;

fn main() {
    let reader = fasta::Reader::from_file("data/test.fa").unwrap();
    for record in reader.records() {
        let record = record.unwrap();
        println!("{}", sequence_to_string(&record.seq()));
    }
}
