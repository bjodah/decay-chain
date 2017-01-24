Decay chain integration using odeint
====================================

.. image:: http://hera.physchem.kth.se:9090/api/badges/bjodah/decay-chain/status.svg
   :target: http://hera.physchem.kth.se:9090/bjodah/decay-chain
   :alt: Build status

A demo where `odeint <http://www.odint.com>`_ is used to integrate a ODE system corrsponding
to a decay chain (cf. radioactive decay without branching).

Example
=======
Run make to run simple tests:

::

   $ make
   g++ -std=c++11 -g -Wall -Wno-unused-local-typedefs -I. -O1  -DVALUE_TYPE_IDX=0 -o decay_chain_dp decay_chain.cpp
   g++ -std=c++11 -g -Wall -Wno-unused-local-typedefs -I. -O1  -DVALUE_TYPE_IDX=1 -o decay_chain_cpp decay_chain.cpp
   g++ -std=c++11 -g -Wall -Wno-unused-local-typedefs -I. -O1  -DVALUE_TYPE_IDX=2 -o decay_chain_mpfr decay_chain.cpp -lmpfr
   g++ -std=c++11 -g -Wall -Wno-unused-local-typedefs -I. -O1  -DVALUE_TYPE_IDX=3 -o decay_chain_gmp decay_chain.cpp -lgmp
   ./decay_chain_dp
   ./decay_chain_cpp
   ./decay_chain_gmp
   # lu decomposition fails with mpfr:
   #./decay_chain_mpfr
   $ echo $?  # zero indicates success
   0
   
for more data you may want to run ``run_all.sh``

::

   $ ./run_all.sh
   + N=27
   + p=1
   + make -B EXTRA_CXXFLAGS=-DVERBOSE decay_chain_dp decay_chain_cpp
   g++ -std=c++11 -g -Wall -Wno-unused-local-typedefs -I. -O1 -DVERBOSE -DVALUE_TYPE_IDX=0 -o decay_chain_dp decay_chain.cpp
   g++ -std=c++11 -g -Wall -Wno-unused-local-typedefs -I. -O1 -DVERBOSE -DVALUE_TYPE_IDX=1 -o decay_chain_cpp decay_chain.cpp
   + method_names=(rosenbrock dopri5 bulirsch_stoer)
   + for method in '{0,1,2}'
   + ./decay_chain_cpp -12 -12 0 -14 1 27 1 27 0
   + bzip2
   log10abstol -12
   log10reltol -12
   log10tend   0
   log10dx     -14
   adaptive    1
   N           27
   p           1
   a           27
   method      0
   
   nrhs 6282
   njac 1047
   4.40088e-14 6.32502e-14 6.35261e-14 4.97635e-14 2.61151e-14 -3.94883e-15 -3.75591e-14 -7.23727e-14 -1.06507e-13 -1.38477e-13 -1.67147e-13 -1.91677e-13 -2.11482e-13 -2.26195e-13 -2.35629e-13 -2.39752e-13 -2.38656e-13 -2.32536e-13 -2.21668e-13 -2.06392e-13 -1.87094e-13 -1.64196e-13 -1.38141e-13 -1.09384e-13 -7.83834e-14 -4.55949e-14 -1.14637e-14 
   + ./decay_chain_dp -12 -12 0 -14 1 27 1 27 0
   + bzip2
   log10abstol -12
   log10reltol -12
   log10tend   0
   log10dx     -14
   adaptive    1
   N           27
   p           1
   a           27
   method      0
   
   nrhs 6282
   njac 1047
   4.40082e-14 6.32506e-14 6.35225e-14 4.97562e-14 2.61137e-14 -3.9465e-15 -3.75628e-14 -7.23917e-14 -1.06514e-13 -1.38486e-13 -1.6716e-13 -1.91694e-13 -2.11508e-13 -2.26211e-13 -2.3564e-13 -2.39765e-13 -2.38663e-13 -2.32543e-13 -2.21675e-13 -2.06396e-13 -1.871e-13 -1.64209e-13 -1.3815e-13 -1.09392e-13 -7.83887e-14 -4.56076e-14 -1.14735e-14 
   + ./instaplot2.py --savefig rosenbrock.png --file-factors 1,-1 rosenbrock_dp.bz2 rosenbrock_cpp.bz2
   0
   /home/bjorn/.local/lib/python3.4/site-packages/matplotlib/backends/backend_gtk3.py:215: Warning: Source ID 7 was not found when attempting to remove it
     GLib.source_remove(self._idle_event_id)
   + for method in '{0,1,2}'
   + ./decay_chain_cpp -12 -12 0 -14 1 27 1 27 1
   + bzip2
   log10abstol -12
   log10reltol -12
   log10tend   0
   log10dx     -14
   adaptive    1
   N           27
   p           1
   a           27
   method      1
   
   nrhs 2360
   njac 0
   3.62407e-14 4.10621e-14 2.45402e-14 -5.19926e-15 -4.1703e-14 -7.9946e-14 -1.16115e-13 -1.47416e-13 -1.71914e-13 -1.88375e-13 -1.96151e-13 -1.9506e-13 -1.85292e-13 -1.67331e-13 -1.41876e-13 -1.09788e-13 -7.20334e-14 -2.96446e-14 1.63176e-14 6.4793e-14 1.14746e-13 1.65183e-13 2.1517e-13 2.63842e-13 3.10411e-13 3.54171e-13 3.94504e-13 
   + ./decay_chain_dp -12 -12 0 -14 1 27 1 27 1
   + bzip2
   log10abstol -12
   log10reltol -12
   log10tend   0
   log10dx     -14
   adaptive    1
   N           27
   p           1
   a           27
   method      1
   
   nrhs 2360
   njac 0
   3.62492e-14 4.10796e-14 2.45624e-14 -5.16427e-15 -4.16594e-14 -7.98962e-14 -1.1606e-13 -1.47358e-13 -1.71852e-13 -1.88311e-13 -1.96083e-13 -1.94986e-13 -1.85222e-13 -1.67257e-13 -1.41802e-13 -1.09711e-13 -7.19563e-14 -2.95649e-14 1.63949e-14 6.487e-14 1.1482e-13 1.65253e-13 2.15243e-13 2.63919e-13 3.10489e-13 3.54244e-13 3.94577e-13 
   + ./instaplot2.py --savefig dopri5.png --file-factors 1,-1 dopri5_dp.bz2 dopri5_cpp.bz2
   0
   /home/bjorn/.local/lib/python3.4/site-packages/matplotlib/backends/backend_gtk3.py:215: Warning: Source ID 7 was not found when attempting to remove it
     GLib.source_remove(self._idle_event_id)
   + for method in '{0,1,2}'
   + ./decay_chain_cpp -12 -12 0 -14 1 27 1 27 2
   + bzip2
   log10abstol -12
   log10reltol -12
   log10tend   0
   log10dx     -14
   adaptive    1
   N           27
   p           1
   a           27
   method      2
   
   nrhs 1403
   njac 0
   8.75757e-18 -1.53875e-16 6.37733e-17 5.22375e-16 8.74583e-16 8.20902e-16 2.18611e-16 -8.98855e-16 -2.35202e-15 -3.86789e-15 -5.13735e-15 -5.86712e-15 -5.82045e-15 -4.84465e-15 -2.88527e-15 1.14688e-17 3.70598e-15 7.97986e-15 1.25552e-14 1.71167e-14 2.13338e-14 2.48827e-14 2.74656e-14 2.88272e-14 2.87676e-14 2.71448e-14 2.40719e-14 
   + ./decay_chain_dp -12 -12 0 -14 1 27 1 27 2
   + bzip2
   log10abstol -12
   log10reltol -12
   log10tend   0
   log10dx     -14
   adaptive    1
   N           27
   p           1
   a           27
   method      2
   
   nrhs 1443
   njac 0
   1.31622e-16 -5.34295e-16 -8.17488e-16 -3.5822e-16 7.32921e-16 2.12677e-15 3.44863e-15 4.34895e-15 4.5762e-15 3.97772e-15 2.52229e-15 2.75821e-16 -2.60209e-15 -5.90326e-15 -9.39873e-15 -1.28161e-14 -1.58918e-14 -1.84158e-14 -2.01453e-14 -2.09763e-14 -2.07317e-14 -1.94116e-14 -1.69847e-14 -1.34216e-14 -8.62678e-15 -4.38538e-15 2.67321e-15 
   + ./instaplot2.py --savefig bulirsch_stoer.png --file-factors 1,-1 bulirsch_stoer_dp.bz2 bulirsch_stoer_cpp.bz2
   0

note that the plotted errors (generated by instaplot2.py) are not
quite correct since there is an interpolation error too, this will
need to be improved by having ``coupled_decay_dp`` report data at same
time points as those outputed from the adaptive integration in
``coupled_decay_cpp``.
   

License
=======
The source code is Open Source and is released under the very permissive
"simplified (2-clause) BSD license". See ``LICENSE.txt`` for further details.
Contributors are welcome to suggest improvements at https://github.com/bjodah/coupled-decay

Author
======
Björn I. Dahlgren, contact:

- gmail address: bjodah
- kth.se address: bda
