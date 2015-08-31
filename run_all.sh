#!/bin/bash -ex
N=${1:-27} # Number of components
p=${2:-1} # decay constant offset
make
rm *.bz2
method_names=(rosenbrock dopri5 bulirsch_stoer)
for method in {0,1,2}; do
    ./cda_mp -12 -12 0 -14 1 $N $p $N $method | bzip2 >${method_names[$method]}_mp.bz2
    ./cda_dp -12 -12 0 -14 1 $N $p $N $method | bzip2 >${method_names[$method]}_dp.bz2
    ./instaplot2.py --savefig ${method_names[$method]}.png --file-factors 1,-1 \
        ${method_names[$method]}_dp.bz2 ${method_names[$method]}_mp.bz2
done
