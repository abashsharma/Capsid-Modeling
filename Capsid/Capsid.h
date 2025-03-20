#pragma once

#include <numbers>

#include "Utils.h"

namespace capsid
{

using namespace std::numbers;

struct Harmonics
{
    Harmonics() = delete;

    // Construct a new Harmonics object
    Harmonics(u_t_ qpoints, u_t_ nmodes, bool randa = false);

    //tex:$$(N+1)^2 + 1$$
    u_t_ NMode;

    //tex:$$n_t = n_p$$
    u_t_ quadpoints;

    double C0;

    double Ct;

    //tex:$$\vec{a}=\left[A_0^0, A_1^0, A_1^1, B_1^1,...,A_N^0,...,A_N^N,B_N^1,...,B_N^N\right]^T$$
    Vec a;

    Vec radius;

    Vec surface_differentials;

    Vec meancurve;

};

/*!
 * \return Sum over the gaussian curvature.
 */
double Calculate_MeanCurve(Harmonics& harmonics);

/*!
 * Get the step size for the integral
 */
#define Q_STEP(x) std::numbers::pi / (x)

/*!
 * Calculate n and m coefficients for spherical harmonics.
 */
inline void Calculate_NM(auto index, auto& n, auto& m) noexcept
{
    n = static_cast<i_t_>(std::floor(std::sqrt(index)));
    m = (index <= n + static_cast<u_t_>(n) * n) ? index - n * n : n * n + n - index;
}

inline auto EnergyModes(u_t_ harmonics) noexcept
{
    return capsid::saturate_cast<u_t_>(std::pow(harmonics + 1, 2));
}

 // Generate an angle from [ANG_ERROR, pi-ANG_ERROR]
inline auto theta_ang(auto step, auto idx) noexcept
{
    return std::min(std::max(step * idx, ANG_ERROR), pi - ANG_ERROR);
}

// Generate an angle from [ANG_ERROR, 2pi-ANG_ERROR]
inline auto phi_ang(auto step, auto idx) noexcept
{
    //return std::min(std::max(2.0 * step * idx, ANG_ERROR), 2 * pi - ANG_ERROR);
    return 2.0 * step * idx;
}

// Normalize the integral by multiplying it by the total sum(dtheta*dphi)
inline double NormalizeIntegral(double integral, u_t_ qpoints) noexcept
{
    return integral * 2 * pi * pi / (qpoints * qpoints);
}

}
