#include <vector>
#include <functional>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/odeint.hpp>

namespace odeintmp {

    enum class StepType : int { bulirsch_stoer, rosenbrock4, dopri5 };

    using std::placeholders::_1;
    using std::placeholders::_2;
    using std::placeholders::_3;

    using boost::numeric::odeint::integrate_const;
    using boost::numeric::odeint::integrate_adaptive;
    using boost::numeric::odeint::make_dense_output;
    using boost::numeric::odeint::make_controlled;
    using boost::numeric::odeint::rosenbrock4;
    using boost::numeric::odeint::runge_kutta_dopri5;
    using boost::numeric::odeint::bulirsch_stoer;
    using boost::numeric::odeint::bulirsch_stoer_dense_out;

    template <typename value_type>
    struct MPODESys {
        typedef boost::numeric::ublas::vector<value_type> vector_type;
        typedef boost::numeric::ublas::matrix<value_type> matrix_type;
        using sys_type = std::pair<std::function<void(const vector_type &, vector_type &, value_type)>,
                                   std::function<void(const vector_type &, matrix_type &,
                                                      const value_type &, vector_type &)> >;
        std::function<void(const vector_type &, value_type)> obs_cb;
        std::function<void(const vector_type &, vector_type &, const value_type)> rhs_cb;
        const int m_ny;
        sys_type m_sys;

        int m_nrhs, m_njac;
        std::vector<value_type> xout;
        std::vector<value_type> yout;
        MPODESys(int ny, sys_type sys) : m_ny(ny), m_sys(sys) {
            obs_cb = std::bind(&MPODESys<value_type>::obs, this, _1, _2);
            rhs_cb = std::bind(&MPODESys<value_type>::rhs, this, _1, _2, _3);
        }
        void obs(const vector_type &yarr, value_type xval) {
            xout.push_back(xval);
            for(size_t i=0 ; i<yarr.size() ; ++i)
                yout.push_back(yarr[i]);
        }

        virtual void rhs(const vector_type &y, vector_type &dydx, const value_type x) = 0;
        virtual void jac(const vector_type & y, matrix_type &J, const value_type & x, vector_type &dfdx) = 0;

        size_t rosenbrock_const(vector_type y0, value_type x0, value_type xend, value_type dx0){
            m_nrhs = 0; m_njac = 0;
            rosenbrock4<value_type> stepper;
            return integrate_const(stepper, this->m_sys, y0, x0, xend, dx0, obs_cb);
        }
        size_t rosenbrock_adaptive(vector_type y0, value_type x0, value_type xend,
                                   value_type dx0, value_type atol, value_type rtol){
            m_nrhs = 0; m_njac = 0;
            auto stepper = make_dense_output<rosenbrock4<value_type> >(atol, rtol);
            return integrate_adaptive(stepper, this->m_sys, y0, x0, xend, dx0, obs_cb);
        }
        size_t rosenbrock_adaptive_nondense(vector_type y0, value_type x0, value_type xend,
                                            value_type dx0, value_type atol, value_type rtol){
            m_nrhs = 0; m_njac = 0;
            auto stepper = make_controlled<rosenbrock4<value_type> >(atol, rtol);
            return integrate_adaptive(stepper, this->m_sys, y0, x0, xend, dx0, obs_cb);
        }
        size_t dopri5_const(vector_type y0, value_type x0, value_type xend, value_type dx0){
            m_nrhs = 0; m_njac = 0;
            runge_kutta_dopri5< vector_type, value_type > stepper;
            return integrate_const(stepper, rhs_cb, y0, x0, xend, dx0, obs_cb);
        }
        size_t dopri5_adaptive(vector_type y0, value_type x0, value_type xend,
                               value_type dx0, value_type atol, value_type rtol){
            m_nrhs = 0; m_njac = 0;
            auto stepper = make_dense_output<runge_kutta_dopri5< vector_type, value_type >>(atol, rtol);
            return integrate_adaptive(stepper, rhs_cb, y0, x0, xend, dx0, obs_cb);
        }
        size_t dopri5_adaptive_nondense(vector_type y0, value_type x0, value_type xend,
                                        value_type dx0, value_type atol, value_type rtol){
            m_nrhs = 0; m_njac = 0;
            auto stepper = make_controlled<runge_kutta_dopri5< vector_type, value_type > >(atol, rtol);
            return integrate_adaptive(stepper, rhs_cb, y0, x0, xend, dx0, obs_cb);
        }
        size_t bulirsch_stoer_const(vector_type y0, value_type x0, value_type xend, value_type dx0){
            m_nrhs = 0; m_njac = 0;
            bulirsch_stoer_dense_out<vector_type, value_type> stepper;
            return integrate_const(stepper, rhs_cb, y0, x0, xend, dx0, obs_cb);
        }
        size_t bulirsch_stoer_adaptive(vector_type y0, value_type x0, value_type xend,
                                       value_type dx0, value_type atol, value_type rtol){
            m_nrhs = 0; m_njac = 0;
            auto stepper = bulirsch_stoer_dense_out<vector_type, value_type>(atol, rtol);
            return integrate_adaptive(stepper, rhs_cb, y0, x0, xend, dx0, obs_cb);
        }
        size_t bulirsch_stoer_adaptive_nondense(vector_type y0, value_type x0, value_type xend,
                                                value_type dx0, value_type atol, value_type rtol){
            m_nrhs = 0; m_njac = 0;
            auto stepper = bulirsch_stoer<vector_type, value_type>(atol, rtol);
            return integrate_adaptive(stepper, rhs_cb, y0, x0, xend, dx0, obs_cb);
        }

    };

}
