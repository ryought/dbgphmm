use bio::io::fasta;

pub fn test() {
    println!("loading");
    let filename = "test.fa";
    // parse_kmers_and_copy_nums(filename);
}

/*
pub fn parse_kmers_and_copy_nums(filename, k: usize) {
    // TODO parse copy num in id field
    // and returns (kmers, copy_nums) that can be used as an input
    // of DbgPHMM::new.
    let mut reader = fasta::Reader::from_file(filename).unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        for window in record.seq().windows(k) {
            println!("seq: {:?}", window);
            let kmer = Kmer::from(window);
        }
    }
}

pub fn read2() {
    let s = b"ATCGATTCGATCGATTCGAT";
    println!("all {:?}", s);
    for w in s.windows(5) {
        println!("window {:?}", w);
    }
}
*/
