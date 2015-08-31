/*
 * coupled_decay.cpp
 *
 * Solving a coupled decay problem
 * with odeint and multiprecision data types
 *
 * Copyright 2015 Björn Dahlgren
 *
 */


#include <array>
#include <functional>
#include <iostream>
#include <stdexcept> // std::logic_error
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/special_functions/binomial.hpp>

// #include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
// #include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>


#if defined(WITH_MULTIPRECISION)
#include <boost/multiprecision/cpp_dec_float.hpp>
typedef boost::multiprecision::cpp_dec_float_50 value_type;
using boost::multiprecision::pow;
using boost::multiprecision::exp;
using boost::multiprecision::log;
#else
#include <cmath>
typedef double value_type;
using std::pow;
using std::exp;
using std::log;
#endif

using namespace std::placeholders;

using boost::numeric::odeint::integrate_const;
using boost::numeric::odeint::integrate_adaptive;
using boost::numeric::odeint::make_dense_output;
using boost::numeric::odeint::rosenbrock4;
using boost::numeric::odeint::runge_kutta_dopri5;
using boost::numeric::odeint::bulirsch_stoer_dense_out;

typedef boost::numeric::ublas::vector<value_type> vector_type;
typedef boost::numeric::ublas::matrix<value_type> matrix_type;

const value_type k0(value_type(17)/value_type(10));
const value_type k1(value_type(23)/value_type(10));
const value_type k2(value_type(29)/value_type(10));

value_type yi1(int i, int p, int a){
    return boost::math::binomial_coefficient<value_type>(p+i, p) * pow(static_cast<value_type>(a), -1 - static_cast<value_type>(p)) *
        pow((a-1)/static_cast<value_type>(a), static_cast<value_type>(i));
}

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

    size_t dopri5_adaptive(vector_type y0, value_type x0, value_type xend,
                        value_type dx0, value_type atol, value_type rtol){
        nrhs = 0; njac = 0;
        auto stepper = make_dense_output<runge_kutta_dopri5< vector_type, value_type >>(atol, rtol);
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

};



size_t run_integration(int log10abstol,
                       int log10reltol,
                       int log10tend,
                       int log10dx,
                       bool adaptive,
                       int N,
                       int p,
                       int a
){
    vector_type y(N);
    value_type atol = pow(value_type(10), value_type(log10abstol));
    value_type rtol = pow(value_type(10), value_type(log10reltol));
    value_type x0 = 0;
    value_type dx0 = pow(value_type(10), value_type(log10dx));
    value_type xend(value_type(pow(value_type(10), value_type(log10tend))));
    vector_type lmbd(N);
    value_type loga = log(static_cast<value_type>(a));
    for (int i=0; i<N; ++i){
        lmbd[i] = (i + p + 1)*loga;
        y[i] = 0;
    }
    y[0] = 1;
    auto cd = CoupledDecay(lmbd);
#if defined(WITH_MULTIPRECISION)
    std::cout.precision(50);
#else
    std::cout.precision(17);
#endif

#if !defined(NDEBUG)
    std::cerr << "log10abstol "<< log10abstol<< '\n'
              << "log10reltol "<< log10reltol<< '\n'
              << "log10tend   "<< log10tend  << '\n'
              << "log10dx     "<< log10dx    << '\n'
              << "adaptive    " << adaptive  << '\n'
              << std::endl;
#endif
    // auto nsteps = cd.rosenbrock_adaptive(y, x0, xend, dx0, atol, rtol);
    // auto nsteps = cd.dopri5_adaptive(y, x0, xend, dx0, atol, rtol);
    auto nsteps = cd.bulirsch_stoer_adaptive(y, x0, xend, dx0, atol, rtol);
#if !defined(NDEBUG)
    std::cerr << "nrhs "<< cd.nrhs << '\n'
              << "njac "<< cd.njac
              << std::endl;
#endif
    for (size_t i=0; i<cd.xout.size(); ++i){
        std::cout << cd.xout[i] << " ";
        for (size_t j=0; j<cd.ny; ++j){
            std::cout << cd.yout[i*cd.ny+j] << " ";
        }
        std::cout << "\n";
    }
#if !defined(NDEBUG)
    for (size_t j=0; j<cd.ny; ++j){
        auto val = cd.yout[nsteps*cd.ny+j];
        auto ref = yi1(j, p, a);
        std::cerr << val - ref << " ";
    }
    std::cerr << "\n";
#endif
    return nsteps;
}


int main(int argc, char **argv)
{
    // see roberts.py for a more user friendly front end.
    if (argc != 9)
        return 1;
    int log10abstol = atoi(argv[1]);
    int log10reltol = atoi(argv[2]);
    int log10tend = atoi(argv[3]);
    int log10dx = atoi(argv[4]);
    int adaptive = atoi(argv[5]);
    int N = atoi(argv[6]);
    int p = atoi(argv[7]);
    int a = atoi(argv[8]);
    std::cerr << "naccept(?): " << run_integration(log10abstol, log10reltol, log10tend,
                                                   log10dx, adaptive == 1, N, p, a) << std::endl;
    return 0;
}
