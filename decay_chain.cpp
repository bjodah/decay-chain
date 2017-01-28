/*
 * decay_chain.cpp
 *
 * Solving a coupled decay problem with arbitrary number of daughters
 * with odeint and (optionally) multiprecision data types
 *
 * Copyright 2015 Björn Dahlgren
 *
 * Open Source, released under the permissive BSD "2-clause license"
 */

#include <iostream>
#include <stdexcept> // std::logic_error
#include "decay_chain.hpp"

// Compile with -DVALUE_TYPE_IDX=# where #: 0=double, 1=cpp_dec_float, 2=mpfr, 3=gmp
#if !defined(VALUE_TYPE_IDX)
#define VALUE_TYPE_IDX 0
#endif

#if VALUE_TYPE_IDX==0
#include <cmath>
  using value_type = double;
  using std::abs;
  using std::pow;
  using std::log;
#elif VALUE_TYPE_IDX==1
  #include <boost/multiprecision/cpp_dec_float.hpp>
  using value_type = boost::multiprecision::cpp_dec_float_50;
  using boost::multiprecision::abs;
  using boost::multiprecision::pow;
  using boost::multiprecision::log;
#elif VALUE_TYPE_IDX==2
  #include <boost/multiprecision/mpfr.hpp>
  using value_type = boost::multiprecision::mpfr_float_50;
  using boost::multiprecision::abs;
  using boost::multiprecision::pow;
  using boost::multiprecision::log;
#elif VALUE_TYPE_IDX==3
  #include <boost/multiprecision/gmp.hpp>
  using value_type = boost::multiprecision::mpf_float_50;
  using boost::multiprecision::abs;
  using boost::multiprecision::pow;
  using boost::multiprecision::log;
#else
  #error "VALUE_TYPE_IDX needs to be either 0, 1, 2 or 3"
#endif

template<typename value_type>
value_type yi1(int i, int p, int a){
    // Analytic solution (given m_lmbd was set correctly)
    return boost::math::binomial_coefficient<value_type>(p+i, p) *
        pow(static_cast<value_type>(a), -1 - static_cast<value_type>(p)) *
        pow((a-1)/static_cast<value_type>(a), static_cast<value_type>(i));
}

using vector_type = typename odeintmp::MPODESys<value_type>::vector_type;
using matrix_type = typename odeintmp::MPODESys<value_type>::matrix_type;

int run_integration(int log10abstol,
                    int log10reltol,
                    int log10tend,
                    int log10dx,
                    int N,
                    int p,
                    int a,
                    int method=0)
{
    vector_type y(N);
    value_type atol = pow(value_type(10), value_type(log10abstol));
    value_type rtol = pow(value_type(10), value_type(log10reltol));
    value_type x0 = 0;
    value_type dx0 = pow(value_type(10), value_type(log10dx));
    value_type xend(value_type(pow(value_type(10), value_type(log10tend))));
    vector_type lmbd(N);
    value_type loga = log(static_cast<value_type>(a));

    // We choose parameters so that we have an easy to evaluate
    // analytic solution:
    for (int i=0; i<N; ++i){
        lmbd[i] = (i + p + 1)*loga;
        y[i] = 0;
    }
    y[0] = 1;
    auto cd = decay_chain::CoupledDecay<value_type>(lmbd);

#if defined(VERBOSE)
    std::cerr << "log10abstol "<< log10abstol<< '\n'
              << "log10reltol "<< log10reltol<< '\n'
              << "log10tend   "<< log10tend  << '\n'
              << "log10dx     "<< log10dx    << '\n'
              << "N           " << N         << '\n'
              << "p           " << p         << '\n'
              << "a           " << a         << '\n'
              << "method      " << method    << '\n'
              << std::endl;
#endif

    size_t nsteps;
    if (method == 0)
        nsteps = cd.rosenbrock_adaptive(y, x0, xend, dx0, atol, rtol);
    else if (method == 1)
        nsteps = cd.dopri5_adaptive(y, x0, xend, dx0, atol, rtol);
    else if (method == 2)
        nsteps = cd.bulirsch_stoer_adaptive(y, x0, xend, dx0, atol, rtol);
    else if (method == 3)
        nsteps = cd.rosenbrock_adaptive_nondense(y, x0, xend, dx0, atol, rtol);
    else if (method == 4)
        nsteps = cd.dopri5_adaptive_nondense(y, x0, xend, dx0, atol, rtol);
    else if (method == 5)
        nsteps = cd.bulirsch_stoer_adaptive_nondense(y, x0, xend, dx0, atol, rtol);
    if (method == 6)
        nsteps = cd.rosenbrock_const(y, x0, xend, dx0);
    else if (method == 7)
        nsteps = cd.dopri5_const(y, x0, xend, dx0);
    else if (method == 8)
        nsteps = cd.bulirsch_stoer_const(y, x0, xend, dx0);
    else
        throw std::logic_error("Unknown method");

#if defined(VERBOSE)
    std::cerr << "nrhs "<< cd.m_nrhs << '\n'
              << "njac "<< cd.m_njac
              << std::endl;
#endif

#if defined(VERBOSE)
    for (size_t i=0; i<cd.xout.size(); ++i){
        std::cout << cd.xout[i] << " ";
        for (int j=0; j<cd.m_ny; ++j){
            std::cout << cd.yout[i*cd.m_ny+j] << " ";
        }
        std::cout << "\n";
    }
#endif

    for (int j=0; j<cd.m_ny; ++j){
        auto val = cd.yout[nsteps*cd.m_ny+j];
        auto ref = yi1<value_type>(j, p, a); // analytic solution
#if defined(VERBOSE)
        std::cerr << val - ref << " ";
#endif
        if (abs(val - ref) > (atol + rtol*ref))
            return 1; // Integration error too big!
    }
#if defined(VERBOSE)
    std::cerr << "\n";
#endif

    return 0;
}


int main(int argc, char **argv)
{
#if VALUE_TYPE_IDX==0
    std::cout.precision(17);
#else
    std::cout.precision(50);
#endif
    int log10abstol, log10reltol, log10tend, log10dx, N, p, a, method;
    if (argc != 9){
        log10abstol = -8;
        log10reltol = -8;
        log10tend = 0;
        log10dx = -14;
        N = 27;
        p = 1;
        a = 27;
        method = 1;
    } else {
        log10abstol = atoi(argv[1]);
        log10reltol = atoi(argv[2]);
        log10tend = atoi(argv[3]);
        log10dx = atoi(argv[4]);
        N = atoi(argv[5]);
        p = atoi(argv[6]);
        a = atoi(argv[7]);
        method = atoi(argv[8]);
    }
    return run_integration(log10abstol, log10reltol, log10tend,
                           log10dx, N, p, a, method);
}
