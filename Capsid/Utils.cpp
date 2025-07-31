#include <random>
#include "Capsid.h"
#include "Utils.h"


namespace capsid
{

std::mt19937& Generator()
{
    thread_local std::random_device rd{};
    thread_local std::mt19937 g(rd());
    return g;
}

Vec randn(u_t_ size)
{
    thread_local std::normal_distribution dist{ 0.0 };
    Vec ret(size);
    for (auto& r : ret)
    {
        r = dist(Generator());
    }
    return ret;
}

Vec randexp(u_t_ size, double a0)
{
    Vec ret{ randn(size) };
    ret[0] = a0;
    double lambda = 3.0;
    for (size_t i = 0; i < size; ++i)
    {
        auto x = static_cast<double>(i) / (size - 1) * a0; // Rescale to [0, a0]
        auto env = std::exp(-lambda * x);
        ret[i] = std::abs(ret[i]) * env;
    }
    return ret;
}

}
