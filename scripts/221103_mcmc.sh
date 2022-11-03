for k in 8 16 24 32 40 48
do
  echo $k
  ./target/release/mcmc -k $k -c 20 -p 0.003 > data221103/tandem_k${k}_c20_p0003.txt
  ./target/release/mcmc -k $k -c 20 -p 0.003 --simple > data221103/simple_k${k}_c20_p0003.txt
done
