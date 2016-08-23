#!/bin/bash

poles_N=( 50 30 20 64 26 22 10 28 10 12 10 6 12 30 20 22 20 11 11 15 12 6 4 60);
ports_Nc=( 3136 6889 60025 1156 2304 8464 1600 2116 1156 27889 1681 25600 1936 576 1225 484 484 784 625 529 1600 19600 5184 10201 );
samples_Ns=( 1001 1228 197 570 690 71 1001 71 1000 27 572 1010 1010 10010 1450 2010 2010 5010 572 572 400 100 110 999 );

:>stats_cpp_all_split_parallel.txt

echo "N;Nc;Ns;iters;RMS_err;time_consumed" > stats_cpp_all_split_parallel.txt

for i in {0..23}
do
    ./vf3 ${poles_N[$i]} ${ports_Nc[$i]} ${samples_Ns[$i]}
done
