
# ./target/release/dbgphmm -k 8 optimize data/r1.fa --true-dbg-fa data/d1.fa -M 10 -V 10 -I 100 -R 0.8 --parallel --start-from-true-copy-nums > data/r1.tsv
# ./target/release/dbgphmm -k 8 stat data/r1.fa > data/r1.json


function from_zero () {
  # initial temperature
  T=100
  # number of iteration
  I=500
  # cooling rate
  R=0.95
  ./target/release/dbgphmm -k 8 optimize data/r1.fa --true-dbg-fa data/d1.fa -M 10 -V 10 -T $T -I $I -R $R --parallel > data/r1_f0_T${T}_R${R}_I${I}.tsv
}

function from_true () {
  # initial temperature
  T=1000
  # number of iteration
  I=500
  # cooling rate
  R=1.0
  # ./target/release/dbgphmm -k 8 optimize data/r1.fa --true-dbg-fa data/d1.fa -M 10 -V 10 -T $T -I $I -R $R --parallel --start-from-true-copy-nums > data/r1_ft_T${T}_R${R}_I${I}.tsv
  python scripts/plotter.py data/r1_ft_T${T}_R${R}_I${I}.tsv data/r1.json
}

from_true
