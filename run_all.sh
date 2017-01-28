#!/bin/bash -ex
N=${1:-27} # Number of components
p=${2:-1} # decay constant offset
make -B EXTRA_CXXFLAGS=-DVERBOSE decay_chain_dp decay_chain_cpp decay_chain_gmp decay_chain_mpfr analytic
method_names=(rosenbrock dopri5 bulirsch_stoer)
make -B analytic
for method in {0,1,2}; do
    for impl in mpfr gmp cpp dp; do
        ./decay_chain_$impl -12 -12 0 -14 $N $p $N $method >${method_names[$method]}_$impl.txt
        ./analytic $N $p $N 1 < ${method_names[$method]}_$impl.txt >${method_names[$method]}_${impl}_err.txt
        ./instaplot2.py --savefig ${method_names[$method]}_${impl}_err.png ${method_names[$method]}_${impl}_err.txt
        ./analytic $N $p $N 0 < ${method_names[$method]}_$impl.txt >${method_names[$method]}_${impl}_ref.txt
        ./instaplot2.py --savefig ${method_names[$method]}_$impl.png --file-factors 1,-1 \
                        ${method_names[$method]}_$impl.txt ${method_names[$method]}_${impl}_ref.txt
    done

    for impl in mpfr gmp cpp; do
        ./decay_chain_$impl -14 -14 0 -30 $N $p $N $method >${method_names[$method]}_${impl}_14.txt
        ./analytic $N $p $N 1 < ${method_names[$method]}_${impl}_14.txt >${method_names[$method]}_${impl}_err_14.txt
        ./instaplot2.py --savefig ${method_names[$method]}_${impl}_err_14.png ${method_names[$method]}_${impl}_err_14.txt
    done

done
./instaplot2.py --savefig analytic.png dopri5_cpp_ref.txt
