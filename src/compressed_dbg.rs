use crate::cycles;
use crate::dbg::{DbgHash, DBG};
use crate::kmer::kmer::{null_kmer, Kmer};
use fnv::FnvHashMap as HashMap;

/// indexed and compressed dbg
/// O(1) access to childs, parents, emissions of index
/// back to kmer content if necessary
/// and it stores indexed cycles
pub struct CompressedDBG {
    k: usize,
    n_kmers: usize,
    childs: Vec<Vec<usize>>,
    parents: Vec<Vec<usize>>,
    kmers: Vec<Kmer>,
    ids: HashMap<Kmer, usize>,
    emissions: Vec<u8>,
    cycles: Vec<Vec<usize>>,
}

impl CompressedDBG {
    /// Construct CompressedDBG from DbgHash
    /// - assign an index to each kmers and convert hashmap into vector
    /// - calc cycles and store in vector
    pub fn from(dbg: &DbgHash, k: usize) -> CompressedDBG {
        // assign an index to each kmers
        let kmers: Vec<Kmer> = dbg.kmers();
        let n_kmers = kmers.len();
        let mut ids: HashMap<Kmer, usize> = HashMap::default();
        for (i, kmer) in kmers.iter().enumerate() {
            ids.insert(kmer.clone(), i);
        }

        let mut childs: Vec<Vec<usize>> = Vec::new();
        let mut parents: Vec<Vec<usize>> = Vec::new();
        let mut emissions: Vec<u8> = Vec::new();
        for kmer in kmers.iter() {
            childs.push(
                dbg.childs(kmer)
                    .iter()
                    .map(|kmer| *ids.get(kmer).unwrap())
                    .collect(),
            );
            parents.push(
                dbg.parents(kmer)
                    .iter()
                    .map(|kmer| *ids.get(kmer).unwrap())
                    .collect(),
            );
            emissions.push(kmer.last());
        }

        let root = null_kmer(k - 1);
        let s = cycles::DbgTree::new(dbg, &root);
        let cycles: Vec<Vec<usize>> = s
            .cycle_keys()
            .iter()
            .map(|key| {
                s.cycle_components(key)
                    .iter()
                    .map(|kmer| *ids.get(kmer).unwrap())
                    .collect()
            })
            .collect();

        CompressedDBG {
            k,
            n_kmers,
            childs,
            parents,
            kmers,
            ids,
            emissions,
            cycles,
        }
    }
    // TODO add construction method from k-1mer compresseddbg
    pub fn n_kmers(&self) -> usize {
        self.n_kmers
    }
    pub fn childs(&self, v: usize) -> &[usize] {
        &self.childs[v]
    }
    pub fn parents(&self, v: usize) -> &[usize] {
        &self.parents[v]
    }
    pub fn is_adjacent(&self, v: usize, w: usize) -> bool {
        self.childs(v).iter().any(|&u| u == w)
    }
    pub fn emission(&self, v: usize) -> u8 {
        self.emissions[v]
    }
    pub fn kmer(&self, v: usize) -> &Kmer {
        &self.kmers[v]
    }
}

pub fn test() {
    let mut d = DbgHash::new();
    // ATCGATTCGAT
    d.add_seq(b"ATCGATTCGATCGATTCGATAGATCG", 8);
    // d.add_seq(b"ATCTCGATCTTGATAGATCG", 8);
    // println!("{}", d.as_dot());
    let cdbg = CompressedDBG::from(&d, 8);
    println!("#cycles: {}", cdbg.cycles.len());
    println!("cycle: {:?}", cdbg.cycles[0]);
}
