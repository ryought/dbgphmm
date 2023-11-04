# dbgphmm

Bayesian genome assembler using de Bruijn graph and profile HMM

## Build

Requires nightly rust

```
cargo build --release
./target/release/dbgphmm -h
```

## Infer

To infer the L-DBG of the read length L, first we construct the initial k0-DBG from the reads by the `draft` subcommand.

```
% ./target/release/dbgphmm draft --help
dbgphmm-draft
Construct draft DBG from reads

USAGE:
    dbgphmm draft [OPTIONS] -k <K> -M <MIN_DEADEND_COUNT> -G <GENOME_SIZE> --dbg-output <DBG_OUTPUT> <READ_FASTA>

ARGS:
    <READ_FASTA>    Input read FASTA filename

OPTIONS:
    -d, --dbg-output <DBG_OUTPUT>    Output dbg filename
    -g, --gfa-output <GFA_OUTPUT>    Output GFA of dbg filename
    -G <GENOME_SIZE>                 Expected size (total number of bases) of target genome
    -h, --help                       Print help information
    -k <K>                           k of DBG
    -m <MIN_COUNT>                   Minimum occurrence of k-mers in read [default: 2]
    -M <MIN_DEADEND_COUNT>           Minimum occurrence of deadend k-mers in read
    -p <P_ERROR>                     Expected error rate of reads. If not specified, it will use
                                     p=0.001 (0.1%) HiFi error rate [default: 0.001]
    -P <N_HAPLOTYPES>                Expected number of haplotypes in target genome if known
```

Then infer the k-DBG (k=k0,...,L) by `infer` subcommand.

```
% ./target/release/dbgphmm infer --help
dbgphmm-infer

USAGE:
    dbgphmm infer [OPTIONS] --dbg-input <DBG_INPUT> --output-prefix <OUTPUT_PREFIX> -K <K_MAX> -G <GENOME_SIZE> -S <GENOME_SIZE_SIGMA> <READ_FASTA>

ARGS:
    <READ_FASTA>    Input read FASTA filename

OPTIONS:
    -c <MAX_CYCLE_SIZE>
            [default: 1000]

    -d, --dbg-input <DBG_INPUT>
            Filename of initial DBG

    -e <P_INFER>
            Error rate of reads while inference [default: 0.00001]

    -G <GENOME_SIZE>
            Expected size (total number of bases) of target genome

    -h, --help
            Print help information

    -I <MAX_ITER>
            Maximum number of iteration of posterior sampling of single k [default: 50]

    -K <K_MAX>
            Target k of DBG

    -o, --output-prefix <OUTPUT_PREFIX>
            Prefix of output files used as a working directory

    -p <P_ERROR>
            Expected error rate of reads. If not specified, it will use p=0.001 (0.1%) HiFi error
            rate [default: 0.001]

        --p0 <P0>
            If probability of copy number being zero is above p0, the k-mer will be regarded as
            zero-copy and purged. In default, 0.8 will be used. Larger p0, bigger the graph
            [default: 0.8]

    -S <GENOME_SIZE_SIGMA>
            Expected size (total number of bases) sigma (standard deviation) of target genome
```

To infer DBG with HiFi reads with the expected genome size `G` and the std var of the genome size `S`, we recommend the following settings.

```
export OMP_NUM_THREADS=1
./target/release/dbgphmm draft -k 40 -G $G -m 2 -M 5 -p 0.0003 -P 2 -d out.dbg reads.fa
./target/release/dbgphmm infer -K 20000 -G $G -S $S -e 0.0003 --p0 0.99 -d out.dbg -o out reads.fa
```

## Simulation experiments

```
# generate synthetic genome and reads
./target/release/draft -h
# infer and evaluate by the true genome
./target/release/infer -h
```

To see the actual usage, see `./scripts/sim.sh`.

## Outputs

### Intermediate results

* `${prefix}.k${k}.gfa`
    GFA output of the inferred DBG and copy numbers
* `${prefix}.k${k}.inspect`
    INSPECT file describing the sampled copy numbers
* `${prefix}.k${k}.dbg`
    serialized DBG
* `${prefix}.k${k}.post`
    (for internal use) dump of posterior distribution

### Final
* `${prefix}.final.gfa`
    Final L-DBG.
    This corresponds to unitig graphs.
* `${prefix}.final.euler.fa`
    A candidate genome corresponding to an Eulerian circuit.
    This corresponds to pseudo-haplotype contigs.


## File formats

### INSPECT

```
# comment
# G section: graph property summary
10000   G       n_edges_full    103992
10000   G       n_edges_compact 7
10000   G       n_nodes_full    103990
10000   G       n_nodes_compact 5
10000   G       n_emittable_edges       86006
10000   G       degree_stats    {(1, 1): 103986, (2, 1): 2, (1, 2): 2}
# C section: sampled copy number assignments and scores
# size of k-mer
#               sample id
#                       P(X|R) posterior probability
#                              P(R|X) likelihood
#                                                P(G) of genome size
#                                                                |G_D(X)| the number of eulerian circuits
#                                                                            |G| genome size
#                                                                                  difference to true copy numbers (if available)
#                                                                                      history of applied update cycles
#                                                                                                   copy number assignment
#                                                                                                          Optional data in JSON
10000   C       0       0.8    -14622.5404367    -9.43613172     0.001       800   0   [S(e2-e5-)]  [2,1]  {}
10000   C       1       0.2    -14622.0570281    -11.4357317     0.002       900   2   []           [2,2]  {}
# ...
# E section: posterior probabilities of each edges
# size of k-mer
#               edge id
#                       true copy number (if available)
#                               expected copy number
#                                       P(X(e)=X*(e)|R) i.e., probability of copy number being true
#                                               P(X(e)=0|R) i.e., probability of copy number being zero
#                                                       marginalized posterior distribution
10000   E       e0      2       2.00000 0.999   0.00000 p(1)=0.000,p(2)=1.000,p(3)=0.000
10000   E       e1      1       1.00000 0.999   0.00000 p(0)=0.000,p(1)=1.000,p(2)=0.000
# ...
```

### DBG

```
# #: comment
# K section: k of DBG
K       4
# N section: node = (k-1)-mer
#       id      sequence of k-1-mer
N       0       nnn
N       1       CAG
# E section: edge = simple path of k-mers
#       id      source  target  joined sequence of k-mers
#                                               current copy number
#                                                       corresponding base ids (internal representation)
E       0       1       0       CAGGAAnnn       1       9,10,11,12,13,14
E       1       1       1       CAGCAG  3       6,7,8
E       2       0       1       nnnTCCCAG       1       0,1,2,3,4,5
```

## Citation

Bayesian genome assembly of segmental duplications by inferring k-mer copy numbers in de Bruijn graphs, WIP


## WIP: python binding

```
source .env/bin/activate
maturin build
pip install --force-reinstall target/wheels/dbgphmm-0.1.0-cp310-cp310-macosx_11_0_arm64.whl
python -c 'import dbgphmm; print(repr(dbgphmm.sum_as_string(1, 2)));'
```
