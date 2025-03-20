#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>
#include <numbers>
#include <limits>
#include <concepts>
#include <stdexcept>
#include <utility>
#include <random>

// Angular error of 1deg
//#define ANG_ERROR (std::numbers::pi / 180.0)
#define ANG_ERROR .001

using u_t_ = uint64_t; // Unsigned integer type. Can probably decrease type size if you want
using i_t_ = int64_t;  // Signed integer type. 


// 1D Vector
using Vec = std::vector<double>;

namespace capsid
{

// Return infinity for a specified type
template<std::floating_point T>
constexpr inline T inf() noexcept
{
    return std::numeric_limits<T>::infinity();
}

// Floating point equality comparison.
template<std::floating_point T>
constexpr inline bool feq(T f1, T f2) noexcept
{
    return std::abs(f1 - f2) <= std::numeric_limits<T>::epsilon();
}

/*!
 * Cosine function that handles pi/2 better.
 */
inline auto cos(auto x) noexcept
{
    const auto c = std::cos(x);
    return capsid::feq(c, 0.0) ? 0.0 : c;
}

// Sine function that handles 0 and n*pi better
inline auto sin(auto x) noexcept
{
    const auto s = std::sin(x);
    return capsid::feq(s, 0.0) ? 0.0 : s;
}

// Compute n!
template<std::integral T>
constexpr inline T factorial(T n) noexcept
{
    // Precompute the first 10 for speed (default N is 9)
    switch (n)
    {
    case 0:
    case 1:
        return 1;
    case 2:
        return 2;
    case 3:
        return 6;
    case 4:
        return 24;
    case 5:
        return 120;
    case 6:
        return 720;
    case 7:
        return 5040;
    case 8:
        return 40320;
    case 9:
        return 362880;
    case 10:
        return 3628800;
    default:
        //tex:$$\Gamma(n) = (n-1)!$$
        return static_cast<T>(std::tgamma(n + 1));
    }
}

// Cross product of two 3D vectors
template<typename T, typename Vector = std::array<T,3>>
inline Vector CrossProduct(const Vector& x, const Vector& y) noexcept
{
    return Vector{
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0]
    };
}

// Compute the normalization constant (length) of two vectors.
template<typename T>
inline double Length(const T& vec) noexcept
{
    return std::sqrt(std::reduce(vec.begin(), vec.end(), 0.0, [](auto sum, auto v) {
        return sum + v * v;
    }));
}

// Compute the dot product of two 3D vectors
template<typename T>
inline double Dot(const T& x, const T& y)
{
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

// Normalize a vector by dividing it by a normalization constant
template<typename T>
inline std::array<T,3> Normalize(std::array<T, 3> vec, T norm)
{
    for (auto& v : vec)
    {
        v /= norm;
    }
    return vec;
}

// Convert degress to radians
template<typename T>
constexpr inline double DegToRad(T deg) noexcept
{
    return deg * std::numbers::pi / 180.0;
}

// Convert radians to degrees
template<typename T>
constexpr inline double RadToDeg(T rad) noexcept
{
    return rad * 180.0 / std::numbers::pi;
}

// std::cmp_less that works on floating point values
template<typename T, typename U>
constexpr bool cmp_less(T t, U u)
{
    if constexpr (std::is_floating_point_v<T> && std::is_floating_point_v<U>)
        if constexpr (sizeof(T) >= sizeof(U))
            return t < static_cast<T>(u);
        else
            return static_cast<U>(t) < u;
    else if constexpr (std::is_floating_point_v<T>)
        return t < static_cast<T>(u);
    else if constexpr (std::is_floating_point_v<U>)
        return static_cast<U>(t) < u;
    else
        return std::cmp_less(t, u);
}

// std::cmp_greater that works on floating point values
template<typename T, typename U>
constexpr bool cmp_greater(T t, U u)
{
    return cmp_less(u, t);
}

// Throws an exception if the value being casted doesn't fit in the new type. 
// static_cast otherwise.
template<typename T, typename U>
inline T saturate_cast(U v)
{
    // Forgo checks in release
#ifndef MAX_SPEED
    constexpr auto tmin = std::numeric_limits<T>::lowest();
    constexpr auto tmax = std::numeric_limits<T>::max();
    constexpr auto umin = std::numeric_limits<U>::lowest();
    constexpr auto umax = std::numeric_limits<U>::max();
    if constexpr (capsid::cmp_less(umin, tmin))
    {
        if (v < U(tmin))
        {
            throw std::underflow_error("Value < U min");
        }
    }
    if constexpr (capsid::cmp_greater(umax, tmax))
    {
        if (v > U(tmax))
        {
            throw std::overflow_error("Value > U max");
        }
    }
#endif
    return static_cast<T>(v);
}

// Reset a vector, e.g. fill it with zeros.
inline void ResetVec(Vec& v) noexcept
{
    std::fill(v.begin(), v.end(), 0);
}

// Get the index for a 1D array as if it's a 2D array.
inline u_t_ Index2D(u_t_ th_idx, u_t_ ph_idx, u_t_ qpoints) noexcept
{
    return th_idx * qpoints + ph_idx;
}

// Return a random number generator engine.
// Thread-safe function.
std::mt19937& Generator();

// Return a vector of n normally distributed random numbers N(0,1)
Vec randn(u_t_ size);

} // namespace capsid

