
# ./target/release/dbgphmm -k 8 optimize data/r1.fa --true-dbg-fa data/d1.fa -M 10 -V 10 -I 100 -R 0.8 --parallel --start-from-true-copy-nums > data/r1.tsv
# ./target/release/dbgphmm -k 8 stat data/r1.fa > data/r1.json


function from_zero () {
  # initial temperature
  T=100
  # number of iteration
  I=1000
  # cooling rate
  R=0.95
  ./target/release/dbgphmm -k 8 benchmark data/r1.fa data/d1.fa -V 10 -T $T -I $I -R $R --parallel annealer > data/r1_f0_T${T}_R${R}_I${I}.tsv
}

function from_true () {
  # initial temperature
  T=1000
  # number of iteration
  I=500
  # cooling rate
  R=1.0
  ./target/release/dbgphmm -k 8 benchmark data/r1.fa data/d1.fa -V 10 -T $T -I $I -R $R --parallel --start-from-true-copy-nums annealer > data/r1_ft_T${T}_R${R}_I${I}.tsv
  python scripts/plotter.py data/r1_ft_T${T}_R${R}_I${I}.tsv data/r1.json
}

function grad () {
  # true
  # ./target/release/dbgphmm -k 8 benchmark data/r1.fa data/d1.fa -V 10 --start-from-true-copy-nums --parallel grad > data/r1.grad.tsv
  ./target/release/dbgphmm -k 8 benchmark data/r1.fa data/d1.fa -V 10 --init-state true --parallel grad -I 10 > data/r1_ft.grad.tsv

  # random
  # ./target/release/dbgphmm -k 8 benchmark data/r1.fa data/d1.fa -V 10 --init-state random --parallel \
  #   grad -I 10 --n-trial 10 --n-basis 3 > data/r1_fr_b3_t50.grad.tsv

  # random
  ./target/release/dbgphmm -k 8 benchmark data/r1.fa data/d1.fa -V 10 --init-state random --parallel grad -I 20 --n-trial 50 --n-basis 10 > data/r1_fr_b10_t50.grad.tsv

  # # read
  # ./target/release/dbgphmm -k 8 benchmark data/r1.fa data/d1.fa -V 1000 --init-state read-count --parallel \
  #   grad -I 10 > data/r1_frc.grad.tsv
}

# from_true
# from_zero
grad
