//#include <array>
#include <functional>
#include <boost/math/special_functions/binomial.hpp>

#include "odeintmp.hpp"

namespace decay_chain{
    using std::placeholders::_1;
    using std::placeholders::_2;
    using std::placeholders::_3;
    using std::placeholders::_4;

    template<typename value_type>
    class CoupledDecay : public odeintmp::MPODESys<value_type> {
    public:
        using vector_type = typename odeintmp::MPODESys<value_type>::vector_type;
        using matrix_type = typename odeintmp::MPODESys<value_type>::matrix_type;

        vector_type m_lmbd;

        CoupledDecay(vector_type lmbd) :
            odeintmp::MPODESys<value_type>(
                lmbd.size(), std::make_pair(std::bind(&CoupledDecay::rhs, this, _1, _2, _3),
                                            std::bind(&CoupledDecay::jac, this, _1, _2, _3, _4))),
            m_lmbd(lmbd) {}
        void rhs(const vector_type &y, vector_type &dydx, const value_type /* x */) override
        {
            for (int i=0; i < (this->m_ny); ++i){
                dydx[i] = -m_lmbd[i]*y[i];
                if (i>0)
                    dydx[i] += m_lmbd[i-1]*y[i-1];
            }
            ++(this->m_nrhs);
        }
        void jac(const vector_type & /* y */, matrix_type &J,
                 const value_type & /* x */, vector_type &dfdx) override {
            for (int i=0; i < (this->m_ny); ++i){
                J(i, i) = -m_lmbd[i];
                if (i>0)
                    J(i, i-1) = m_lmbd[i-1];
                dfdx[i] = 0;
            }
            ++(this->m_njac);
        }
    };
}
