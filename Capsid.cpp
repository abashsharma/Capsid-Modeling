#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <mdspan>
#include <numeric>
#include <array>
#include <cstdlib>
#include <print>
#include <iostream>
#include <numbers>

#include "Defs.h"
#include "Utils.h"
#include "Capsid.h"

namespace capsid
{

using namespace std::numbers;

#define NM_FROM_INDEX(i) \
    i_t_ m, n; \
    Calculate_NM(i, n, m)


/*!
 * Index an mdspan.
 */
#define _I(...) std::array{ __VA_ARGS__ }

namespace
{

CAP_FORCEINLINE auto sh_ang(auto m, auto phi) noexcept
{
    return (m < 0) ? capsid::sin(m * phi) : capsid::cos(m * phi);
}

// Condon-Shortley phase term
//tex:$$(-1)^m$$
CAP_FORCEINLINE auto cs_phase(auto m) noexcept
{
    return (m % 2 == 0) ? 1 : -1;
}

inline auto alegendre(auto n, auto m, auto cos_t)
{
    // stl associated legendre omits the phase term
    return cs_phase(std::abs(m)) * std::assoc_legendre(capsid::saturate_cast<unsigned>(n),
        capsid::saturate_cast<unsigned>(std::abs(m)), cos_t);
}

/*!
 * Calculate normalization factor
 */
inline auto Calculate_SphNorm(auto n, auto m) noexcept
{
    //tex:$$f_n^m = \sqrt{\frac{(2n+1)(n-m)!}{4\pi(n+m)!}}$$
    const auto pref = static_cast<double>(capsid::factorial<u_t_>(n - m)) / capsid::factorial<u_t_>(n + m);
    return std::sqrt((2.0 * n + 1) / (4 * pi) * pref);
}

inline auto Calculate_SurfaceHarmonics(auto n, auto m, auto theta, auto mphi) noexcept
{
    //tex: $$Y_n^m(\theta, \phi) = 
    // \begin{cases}
    // f_n^mP_n^m(\mu)cos(m\phi)), & \text{if } m \geq 0,\\
    // f_n^{|m|}P_n^{|m|}(\mu)sin(|m|\phi) & \text{if } m < 0
    // \end{cases}$$
    return mphi * std::sph_legendre(capsid::saturate_cast<unsigned>(n),
        capsid::saturate_cast<unsigned>(std::abs(m)), theta);
}

inline auto del_legendre(auto n, auto m, auto cos_t, auto sin_t) noexcept
{
    //tex: $$\partial_\theta P^m_n(\mu) = \frac{-1}{sin(\theta)}\left((n+1)cos(\theta)P_n^m(\mu) - (n-m+1)P^m_{n+1}(\mu)\right)$$
    return ((n - m + 1.0) * alegendre(n + 1, m, cos_t)
        - (n + 1.0) * cos_t * alegendre(n, m, cos_t)) / sin_t;
}

inline auto del2_legendre(auto n, auto m, auto cos_t, auto sin_t) noexcept
{
    //tex: $$\partial_\theta^2 P_n^m(\mu) = \left((n+1+(n+1)^2cos^2(\theta))P_n^m(\mu) - 2cos(\theta)(n-m+1)(n+2)P_{n+1}^m(\mu)+(n-m+1)(n-m+2)P_{n+2}^m(\mu)\right)\frac{1}{sin^2(\theta)}$$
    const auto n1 = n + 1.0;
    const auto p1 = (n1 + n1 * n1 * cos_t * cos_t) * alegendre(n, m, cos_t);
    const auto p2 = -2 * cos_t * (n1 - m) * (n + 2) * alegendre(n1, m, cos_t);
    const auto p3 = (n1 - m) * (n - m + 2) * alegendre(n + 2, m, cos_t);
    return (p1 + p2 + p3) / (sin_t * sin_t);
}

inline double _calc_r(i_t_ n, i_t_ m, double theta, double mphi, double a) noexcept
{
    return a * Calculate_SurfaceHarmonics(n, m, theta, mphi);
}

inline double _calc_rphi(i_t_ n, i_t_ m, double theta, double mphi, double a) noexcept
{
    return -m * a * Calculate_SurfaceHarmonics(n, m, theta, mphi);
}

inline double _calc_rphiphi(i_t_ n, i_t_ m, double cos_t, double mphi, double a) noexcept
{
    return -m * m * _calc_r(n, m, cos_t, mphi, a);
}

inline double _calc_rtheta(i_t_ n, i_t_ m, double cos_t, double sin_t, double mphi, double a) noexcept
{
    return a * mphi * Calculate_SphNorm(n, std::abs(m)) * del_legendre(n, m, cos_t, sin_t);
}

inline double _calc_rthetatheta(i_t_ n, i_t_ m, double cos_t, double sin_t, double mphi, double a) noexcept
{
    return a * mphi * Calculate_SphNorm(n, std::abs(m)) * del2_legendre(n, m, cos_t, sin_t);
}

inline double _calc_rthetaphi(i_t_ n, i_t_ m, double cos_t, double sin_t, double mphi, double a) noexcept
{
    return -m * a * mphi * Calculate_SphNorm(n, std::abs(m)) * del_legendre(n, m, cos_t, sin_t);
}

}

Harmonics::Harmonics(u_t_ qpoints, u_t_ nmodes, bool randa)
    : NMode(nmodes)
    , quadpoints(qpoints)
    , C0(0.0)
    , Ct(0.0)
    , a(randa ? randn(nmodes) : Vec(nmodes, 0))
    , radius(qpoints*qpoints, 0)
    , surface_differentials(qpoints*qpoints, 0)
    , meancurve(qpoints*qpoints, 0)
{
}


double Calculate_MeanCurve(Harmonics& h)
{
    const double d_quadpoint = Q_STEP(h.quadpoints);
    double K = 0.0;

    for (u_t_ k = 0; k < h.quadpoints; ++k)
    {
        const auto theta = theta_ang(k, d_quadpoint);
        const auto cos_t = capsid::cos(theta);
        const auto sin_t = capsid::sin(theta);

        for (u_t_ l = 0; l < h.quadpoints; ++l)
        {
            const auto phi   = phi_ang(l, d_quadpoint);
            const auto cos_p = capsid::cos(phi);
            const auto sin_p = capsid::sin(phi);

            // r and its derivatives
            double r             = 0.0;
            double r_phi         = 0.0;
            double r_theta       = 0.0;
            double r_phi_phi     = 0.0;
            double r_theta_theta = 0.0;
            double r_theta_phi   = 0.0;
            for (u_t_ i = 0; i < h.NMode; ++i)
            {
                NM_FROM_INDEX(i);
                const auto a   = h.a[i];
                const auto S   = std::sph_legendre(n, capsid::saturate_cast<unsigned>(std::abs(m)), theta);
                const auto f   = Calculate_SphNorm(n, std::abs(m));
                const auto dS  = f * del_legendre( n, std::abs(m), cos_t, sin_t);
                const auto d2S = f * del2_legendre(n, std::abs(m), cos_t, sin_t);
                const auto cmp = capsid::cos(m * phi);
                const auto smp = capsid::sin(m * phi);
                double r_      = 0.0;
                if (m >= 0) [[likely]]
                {
                    r_           = a * cmp * S;
                    r_phi       += a * smp * -m * S;
                    r_theta     += a * cmp * dS;
                    r_theta_phi += a * smp * -m * dS;
                    if (m != 0 || n != 0) [[likely]]
                    {
                        r_theta_theta += a * cmp * d2S;
                    }
                }
                else
                {
                    r_             = a * smp * S;
                    r_phi         += a * cmp * m * S;
                    r_theta       += a * smp * dS;
                    r_theta_phi   += a * cmp * m * dS;
                    r_theta_theta += a * smp * d2S;
                }
                //tex:$$r_{\phi\phi}=-m^2r$$
                r_phi_phi += -(m * m * r_);
                r         += r_;
            }
            
            const auto idx = Index2D(k, l, h.quadpoints);
            h.radius[idx] = r;

            // We get division by 0 when the radius is 0
            // Gaussian curve + mean curve should be 0 anyways I think?
            // Anyways, this really only happens with a small subset of shapes
            // such as the spherical harmonics for Y(2,1)
            if (r != 0.0) [[likely]]
            {

                //tex:$$\vec{x} = \left[r(\theta,\phi)sin(\theta)cos(\phi), r(\theta,\phi)sin(\theta)sin(\phi), r(\theta,\phi)cos(\theta) \right]$$
                // Compute derivatives of x with respect to theta and phi
                std::array<double, 3> xt, xp;
                xt[0] = cos_p * (r * cos_t + r_theta * sin_t);
                xt[1] = sin_p * (r * cos_t + r_theta * sin_t);
                xt[2] = r_theta * cos_t - r * sin_t;

                xp[0] = sin_t * (-r * sin_p + r_phi * cos_p);
                xp[1] = sin_t * (r * cos_p + r_phi * sin_p);
                xp[2] = r_phi * cos_t;

                const auto xc{ capsid::CrossProduct<double>(xt, xp) };

                //tex:$$R=\left|\vec{x_\theta} \times \vec{x_\phi}\right|$$
                const auto R = capsid::Length(xc);

                // A couple precomputed values so that we have less CPU operations
                // Avoid the use of std::pow so that we don't have to jump to
                // another function address
                const auto rp2 = r * r;
                const auto rst = r * sin_t;
                const auto rtp2 = r_theta * r_theta;
                const auto rpp2 = r_phi * r_phi;


                // I order fundamental form coefficients.
                //tex:$$E=r_\theta^2+r^2$$
                const auto E = rtp2 + rp2;
                //tex:$$F=r_\theta r_\phi$$
                const auto F = r_theta * r_phi;
                //tex:$$G=r_\phi^2 + r^2sin^2(\theta)$$
                const auto G = rpp2 + rst * rst;

                //HUGE TODO:
                //************
                // Check to see whether the - signs I placed in front of e and g
                // are correct.
                // I think that they're the proper form for the convention used,
                // but I'm not absolutely sure. There's some funky things going on with
                // shape operators.

                // II order fundamental form coefficients
                //tex:$$e=R^{-1}(-rsin(\theta)(r^2 +2r_\theta^2 -r r_{\theta\theta}))$$
                const auto e = -(-rst * (rp2 + 2 * rtp2 - r * r_theta_theta)) / R;
                //tex:$$f=R^{-1}(r(-r_\phi(rcos(\theta)+2r_\theta sin(\theta))+rr_{\theta\phi} sin(\theta))$$

                const auto f = (-r * (r_phi * (r * cos_t + 2 * sin_t * r_theta) - rst * r_theta_phi)) / R;
                //const auto f = (r * (-r_phi * (r * cos_t + 2 * r_theta * sin_t) + r_theta_phi * rst)) / R;
                //tex:$$g=R^{-1}(rsin(\theta)(-2r_\phi^2+r(-rsin^2(\theta)+r_{\phi\phi}+r_{\theta}sin(\theta)cos(\theta))))$$
                const auto g = -(-rst * (rst * rst + 2 * rpp2 - r * (r_phi_phi + cos_t * sin_t * r_theta))) / R;
                //const auto g = (rst * (-2 * rpp2 + r * (-rst * sin_t + r_phi_phi + r_theta * cos_t * sin_t))) / R;

                const auto HNumer = (e * G - 2 * f * F + g * E);
                const auto KNumer = (e * g - f * f);
                const auto Denom = (E * G - F * F);

                //tex:$$K=\frac{eg-f^2}{EG-F^2}$$
                const auto Kt = KNumer / Denom;

                //tex:$$dS=r*\sqrt{r_\phi^2+(r_\theta sin(\theta))^2 + (r^2 sin(\theta))^2}d\theta d\phi$$
                //$$\qquad=\sqrt{EG-F^2} d\theta d\phi$$

                const auto dS = std::sqrt(Denom);
                //double dS = r * std::sqrt(r_phi * r_phi + std::pow(r_theta * sin_t, 2) + std::pow(r * sin_t, 2));

                //tex:$$\int KdS = 2\pi \chi(M)$$
                K += Kt * dS;

                //tex:$$H=\frac{eG+gE-2fF}{2\left(EG-F^2\right)}$$
                h.meancurve[idx] = HNumer / 2.0 / Denom;
                h.surface_differentials[idx] = dS;
            }
        }
    }
    return NormalizeIntegral(K, h.quadpoints);
}

}
