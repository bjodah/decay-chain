#include <cmath>
#include <iostream>
#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "bateman.hpp"

using Real_t = boost::multiprecision::cpp_dec_float_100;
using vec_t = std::vector<Real_t>;
Real_t exp_cb(Real_t arg){
    return boost::multiprecision::exp(arg);
}

int main(int argc, char *argv[]){
    int N, p, a;
    if (argc != 4){
        std::cerr << "Expected 3 arguments (N, p, a)\n";
        return 1;
    }
    N = atoi(argv[1]);
    p = atoi(argv[2]);
    a = atoi(argv[3]);
    Real_t loga = boost::multiprecision::log(static_cast<Real_t>(a));
    const auto one = static_cast<Real_t>(1);
    vec_t lmbd;

    for (int i=0; i<N; ++i){
        lmbd.push_back((i + p + one)*loga);
    }
#if !defined(NDEBUG)
    std::cout << "N=" << N << ", p=" << p << ", a=" << a << ", loga=" << loga  << ", one=" << one << '\n';
    for (int i=0; i<N; ++i){
        std::cerr << lmbd[i] << " ";
    }
    std::cerr << '\n';
#endif

    vec_t vals;
    for (std::string line; std::getline(std::cin, line); ){
        if (line.size() == 0)
            break;
        std::string first;
        std::getline(std::istringstream(line), first, ' ');
        vals.push_back((Real_t)first);
    }

   std::cout << std::setprecision(50);
   for (size_t i=0; i<vals.size(); ++i){
       std::cout << vals[i] << " ";
       auto p = bateman::bateman_parent(lmbd, vals[i], exp_cb);
       for (auto v : p)
           std::cout << v << " ";
       std::cout << '\n';
   }
   return 0;
}
