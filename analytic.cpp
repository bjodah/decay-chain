#include <cmath>
#include <iostream>
#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "bateman.hpp"

using Real_t = boost::multiprecision::cpp_dec_float_50;
using vec_t = std::vector<Real_t>;
Real_t exp_cb(Real_t arg){
    return boost::multiprecision::exp(arg);
}

int main(int argc, char *argv[]){
    int N, p, a, diff;
    if (argc != 5){
        std::cerr << "Expected 5 arguments (N, p, a, diff)\n";
        return 1;
    }
    N = atoi(argv[1]);
    p = atoi(argv[2]);
    a = atoi(argv[3]);
    diff = atoi(argv[4]);

    Real_t loga = boost::multiprecision::log(static_cast<Real_t>(a));
    const auto one = static_cast<Real_t>(1);
    vec_t lmbd;

    for (int i=0; i<N; ++i){
        lmbd.push_back((i + p + one)*loga);
    }

    if (diff)
        std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 2);
    else
        std::cout << std::setprecision(50);

    vec_t vals;
    int i = 0;
    for (std::string line; std::getline(std::cin, line); ++i){
        if (line.size() == 0)
            break;
        std::string item;
        auto linestream = std::istringstream(line);
        std::getline(linestream, item, ' ');
        vals.push_back((Real_t)item);

        std::cout << vals[i] << " ";
        auto p = bateman::bateman_parent(lmbd, vals[i], exp_cb);
        if (diff) {
            int j = 0;
            for (; std::getline(linestream, item, ' '); ++j)
                std::cout << p[j] - (Real_t)item << " ";
            if (j != N)
                throw std::logic_error("Unexpected line length");
        } else {
            for (auto v : p)
                std::cout << v << " ";
        }
        std::cout << '\n';
    }
    return 0;
}
