#!/bin/bash -ex
N=${1:-27}
p=${2:-1}
make
rm *.txt.bz2
./cda_dp -12 -12 0 -14 1 $N $p $N | bzip2 >dp.txt.bz2
./cda_mp -12 -12 0 -14 1 $N $p $N | bzip2 >mp.txt.bz2
./instaplot2.py --file-factors 1,-1 dp.txt.bz2 mp.txt.bz2
