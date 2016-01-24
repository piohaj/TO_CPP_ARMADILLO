#!/bin/bash

poles_N=( 50 30 20 64 26 22 10 28 10 12 10 6 12 30 20 22 20 11 11 15 12 6 4 60);
ports_Nc=( 56 83 245 34 48 92 40 46 34 167 41 160 44 24 35 22 22 28 25 23 40 140 72 101 );
samples_Ns=( 10010 12280 1970 5700 6900 710 10010 710 10000 270 5720 1010 3010 10010 1450 2010 2010 5010 5720 5720 400 100 110 10000 );

:>stats_cpp_parallel.txt

echo "N;Nc;Ns;iters;RMS_err;time_consumed" > stats_cpp_parallel.txt

for i in {0..23}
do
    ./vf3 ${poles_N[$i]} ${ports_Nc[$i]} ${samples_Ns[$i]}
done
