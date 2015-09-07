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


#include <array>
#include <functional>
#include <iostream>
#include <stdexcept> // std::logic_error
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/special_functions/binomial.hpp>

// Compile with -DVALUE_TYPE_IDX=# where #: 0=double, 1=cpp_dec_float, 2=mpfr, 3=gmp
#if !defined(VALUE_TYPE_IDX)
#define VALUE_TYPE_IDX 0
#endif

#if VALUE_TYPE_IDX==0
#include <cmath>
  typedef double value_type;
  using std::abs;
  using std::pow;
  using std::log;
#elif VALUE_TYPE_IDX==1
  #include <boost/multiprecision/cpp_dec_float.hpp>
  typedef boost::multiprecision::cpp_dec_float_50 value_type;
  using boost::multiprecision::abs;
  using boost::multiprecision::pow;
  using boost::multiprecision::log;
#elif VALUE_TYPE_IDX==2
  #include <boost/multiprecision/mpfr.hpp>
  typedef boost::multiprecision::mpfr_float_50 value_type;
  using boost::multiprecision::abs;
  using boost::multiprecision::pow;
  using boost::multiprecision::log;
#elif VALUE_TYPE_IDX==3
  #include <boost/multiprecision/gmp.hpp>
  typedef boost::multiprecision::mpf_float_50 value_type;
  using boost::multiprecision::abs;
  using boost::multiprecision::pow;
  using boost::multiprecision::log;
#else
  #error "Unknown value type index"
#endif

using namespace std::placeholders;

using boost::numeric::odeint::integrate_const;
using boost::numeric::odeint::integrate_adaptive;
using boost::numeric::odeint::make_dense_output;
using boost::numeric::odeint::make_controlled;
using boost::numeric::odeint::rosenbrock4;
using boost::numeric::odeint::rosenbrock4_controller;
using boost::numeric::odeint::runge_kutta_dopri5;
using boost::numeric::odeint::controlled_runge_kutta;
using boost::numeric::odeint::bulirsch_stoer;
using boost::numeric::odeint::bulirsch_stoer_dense_out;

typedef boost::numeric::ublas::vector<value_type> vector_type;
typedef boost::numeric::ublas::matrix<value_type> matrix_type;

class CoupledDecay {
public:
    size_t nrhs, njac;
    vector_type m_lmbd;
    const size_t ny;
    std::vector<value_type> xout;
    std::vector<value_type> yout;
private:
    std::pair<std::function<void(const vector_type &, vector_type &, value_type)>,
              std::function<void(const vector_type &, matrix_type &,
                                 const value_type &, vector_type &)> > system;
    void obs(const vector_type &yarr, value_type xval) {
        xout.push_back(xval);
        for(size_t i=0 ; i<yarr.size() ; ++i)
            yout.push_back(yarr[i]);
    }
public:

    CoupledDecay(vector_type lmbd) :
        m_lmbd(lmbd), ny(lmbd.size()),
        system(std::make_pair(std::bind(&CoupledDecay::rhs, this, _1, _2, _3),
                              std::bind(&CoupledDecay::jac, this, _1, _2, _3, _4))) {}
    void rhs(const vector_type &y, vector_type &dydx, value_type x)
    {
        for (size_t i=0; i < (this->ny); ++i){
            dydx[i] = -m_lmbd[i]*y[i];
            if (i>0)
                dydx[i] += m_lmbd[i-1]*y[i-1];
        }
        ++nrhs;
    }
    void jac(const vector_type & y, matrix_type &J, const value_type & x, vector_type &dfdx) {
        for (size_t i=0; i < (this->ny); ++i){
            J(i, i) = -m_lmbd[i];
            if (i>0)
                J(i, i-1) = m_lmbd[i-1];
            dfdx[i] = 0;
        }
        ++njac;
    }
    size_t rosenbrock_adaptive(vector_type y0, value_type x0, value_type xend,
                               value_type dx0, value_type atol, value_type rtol){
        nrhs = 0; njac = 0;
        auto stepper = make_dense_output<rosenbrock4<value_type> >(atol, rtol);
        auto obs_cb = std::bind(&CoupledDecay::obs, this, _1, _2);
        return integrate_adaptive(stepper, this->system, y0, x0, xend, dx0, obs_cb);
    }
    size_t rosenbrock_adaptive_nondense(vector_type y0, value_type x0, value_type xend,
                               value_type dx0, value_type atol, value_type rtol){
        nrhs = 0; njac = 0;
        auto stepper = make_controlled<rosenbrock4<value_type> >(atol, rtol);
        auto obs_cb = std::bind(&CoupledDecay::obs, this, _1, _2);
        return integrate_adaptive(stepper, this->system, y0, x0, xend, dx0, obs_cb);
    }

    size_t dopri5_adaptive(vector_type y0, value_type x0, value_type xend,
                        value_type dx0, value_type atol, value_type rtol){
        nrhs = 0; njac = 0;
        auto stepper = make_dense_output<runge_kutta_dopri5< vector_type, value_type >>(atol, rtol);
        auto obs_cb = std::bind(&CoupledDecay::obs, this, _1, _2);
        return integrate_adaptive(stepper, std::bind(&CoupledDecay::rhs, this, _1, _2, _3),
                                  y0, x0, xend, dx0, obs_cb);
    }
    size_t dopri5_adaptive_nondense(vector_type y0, value_type x0, value_type xend,
                        value_type dx0, value_type atol, value_type rtol){
        nrhs = 0; njac = 0;
        auto stepper = make_controlled<runge_kutta_dopri5< vector_type, value_type > >(atol, rtol);
        auto obs_cb = std::bind(&CoupledDecay::obs, this, _1, _2);
        return integrate_adaptive(stepper, std::bind(&CoupledDecay::rhs, this, _1, _2, _3),
                                  y0, x0, xend, dx0, obs_cb);
    }
    size_t bulirsch_stoer_adaptive(vector_type y0, value_type x0, value_type xend,
                        value_type dx0, value_type atol, value_type rtol){
        nrhs = 0; njac = 0;
        auto stepper = bulirsch_stoer_dense_out<vector_type, value_type>(atol, rtol);
        auto obs_cb = std::bind(&CoupledDecay::obs, this, _1, _2);
        return integrate_adaptive(stepper, std::bind(&CoupledDecay::rhs, this, _1, _2, _3),
                                  y0, x0, xend, dx0, obs_cb);
    }
    size_t bulirsch_stoer_adaptive_nondense(vector_type y0, value_type x0, value_type xend,
                              value_type dx0, value_type atol, value_type rtol){
        nrhs = 0; njac = 0;
        auto stepper = bulirsch_stoer<vector_type, value_type>(atol, rtol);
        auto obs_cb = std::bind(&CoupledDecay::obs, this, _1, _2);
        return integrate_adaptive(stepper, std::bind(&CoupledDecay::rhs, this, _1, _2, _3),
                                  y0, x0, xend, dx0, obs_cb);
    }
    static value_type yi1(int i, int p, int a){
        // Analytic solution (given m_lmbd was set correctly)
        return boost::math::binomial_coefficient<value_type>(p+i, p) *
            pow(static_cast<value_type>(a), -1 - static_cast<value_type>(p)) *
            pow((a-1)/static_cast<value_type>(a), static_cast<value_type>(i));
}


};


int run_integration(int log10abstol,
                    int log10reltol,
                    int log10tend,
                    int log10dx,
                    bool adaptive,
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
    auto cd = CoupledDecay(lmbd);

#if defined(VERBOSE)
    std::cerr << "log10abstol "<< log10abstol<< '\n'
              << "log10reltol "<< log10reltol<< '\n'
              << "log10tend   "<< log10tend  << '\n'
              << "log10dx     "<< log10dx    << '\n'
              << "adaptive    " << adaptive  << '\n'
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
    else
        throw std::logic_error("Unknown method");

#if defined(VERBOSE)
    std::cerr << "nrhs "<< cd.nrhs << '\n'
              << "njac "<< cd.njac
              << std::endl;
#endif

#if defined(VERBOSE)
    for (size_t i=0; i<cd.xout.size(); ++i){
        std::cout << cd.xout[i] << " ";
        for (size_t j=0; j<cd.ny; ++j){
            std::cout << cd.yout[i*cd.ny+j] << " ";
        }
        std::cout << "\n";
    }
#endif

    for (size_t j=0; j<cd.ny; ++j){
        auto val = cd.yout[nsteps*cd.ny+j];
        auto ref = CoupledDecay::yi1(j, p, a); // analytic solution
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
    int log10abstol, log10reltol, log10tend, log10dx, adaptive, N, p, a, method;
    if (argc != 10){
        log10abstol = -12;
        log10reltol = -12;
        log10tend = 0;
        log10dx = -14;
        adaptive = 1;
        N = 27;
        p = 1;
        a = 27;
        method = 0;
    } else {
        log10abstol = atoi(argv[1]);
        log10reltol = atoi(argv[2]);
        log10tend = atoi(argv[3]);
        log10dx = atoi(argv[4]);
        adaptive = atoi(argv[5]);
        N = atoi(argv[6]);
        p = atoi(argv[7]);
        a = atoi(argv[8]);
        method = atoi(argv[9]);
    }
    return run_integration(log10abstol, log10reltol, log10tend,
                           log10dx, adaptive == 1, N, p, a, method);
}
