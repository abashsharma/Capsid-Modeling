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

}
