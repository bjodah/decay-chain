#!/bin/bash -ex
N=${1:-27} # Number of components
p=${2:-1} # decay constant offset
make -B EXTRA_CXXFLAGS=-DVERBOSE decay_chain_dp decay_chain_cpp decay_chain_gmp decay_chain_mpfr analytic
method_names=(rosenbrock dopri5 bulirsch_stoer)
make -B analytic
for method in {0,1,2}; do
    for impl in mpfr gmp cpp dp; do
        ./decay_chain_$impl -12 -12 0 -14 $N $p $N $method >${method_names[$method]}_$impl.txt
        ./analytic $N $p $N ${method_names[$method]}_$impl.txt >${method_names[$method]}_ref.txt
        ./instaplot2.py --savefig ${method_names[$method]}.png --file-factors 1,-1 \
                        ${method_names[$method]}_$impl.txt ${method_names[$method]}_ref.txt
    done
done
