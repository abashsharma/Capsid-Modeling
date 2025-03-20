#include <stdexcept>
#include <numbers>
#include <string>
#include <random>
#include <utility>
#include <algorithm>
#include <iostream>
#include "CSV.h"
#include "Files.h"
#include "Utils.h"
#include "Capsid.h"
#include "Optimize.h"

namespace
{

//Todo: Implement ncg
void ncg(capsid::Harmonics& h)
{

}

void monte_carlo(capsid::Harmonics& h)
{
    //Todo: Remove CSV file writing
    CSVWriter Ewrite("e.csv");

    const auto fname{ "mc_capsid_" + std::to_string(h.C0) + "_" + std::to_string(h.Ct) + ".xyz" };

    constexpr auto NSAMPLES = 1000;
    const auto kT = .0005; 

    std::uniform_real_distribution<double> dAdist(-.01, .01);

    std::uniform_int_distribution<> accept(0, 100);

    std::uniform_int_distribution<> toPerturb(1, h.a.size());

    auto Kcalc = Calculate_MeanCurve(h);
    std::cout << "K: " << Kcalc << '\n';
    auto E = std_bending_energy(h, Kcalc);
    capsid::SaveRadii(h, fname);
    Vec a0{ std::move(h.a) };

    Ewrite.write(E);

    // How 2 converge???
    for (int i = 0; i < NSAMPLES; ++i)
    {
        int ntb = toPerturb(capsid::Generator());
        Vec anew{ a0 };
        for (int j = 0; j < ntb; ++j)
        {
            auto idx   = toPerturb(capsid::Generator()) - 1;
            anew[idx] += dAdist(capsid::Generator());
        }
        // set a subset of as from a distribution
        h.a = std::move(anew);
        Kcalc  = Calculate_MeanCurve(h);
        const double Enew = std_bending_energy(h, Kcalc);

        const auto dE = Enew - E;
        // match kT to tolerance
        // tune kt based on typical energy variance as make proposal
        // set kt within 1 stdev of a typical dE
        if (dE < 0  || accept(capsid::Generator()) / 100.0 < std::exp(-dE / kT))
        {
            capsid::SaveRadii(h, fname, true);
            a0 = std::move(h.a);
            E  = Enew;
            Ewrite.write(E);
        }
    }
    h.a = std::move(a0);
    // Recalculate radii
    capsid::Calculate_MeanCurve(h);
    capsid::SaveRadii(h, fname, true);
}

}

void optimize(capsid::Harmonics& h, std::string_view method)
{
    if (method.compare("NCG") == 0)
    {
        ncg(h);
    }
    else if (method.compare("MC") == 0)
    {
        monte_carlo(h);
    }
    else
    {
        throw std::invalid_argument("Unrecognized optimization method");
    }
}

//tex:$$E\left[\Gamma(\vec{x})\right] = \int\left(\frac{1}{2}\mathbb{K}_C(2H-C_0)^2+\mathbb{K}_GK\right)dS$$
//$$\qquad\qquad=\frac{1}{2}\mathbb{K}_C\int(2H-C_0)^2dS + \mathbb{K}_G\int KdS$$
double std_bending_energy(const capsid::Harmonics& h, double K)
{
    //Todo: Set KC and KG
    constexpr auto KC = 1.0; // Bending modulus
    // The authors set KG to 0 in their code???
    constexpr auto KG = 1.0; // Gaussian saddle-splay modulus
    const auto size   = h.quadpoints * h.quadpoints;
    double e_sum      = 0.0;
    for (u_t_ i = 0; i < size; ++i)
    {
        const auto diff = 2 * h.meancurve[i] - h.C0;
        e_sum += diff * diff * h.surface_differentials[i];
    }

    return .5 * KC * capsid::NormalizeIntegral(e_sum, h.quadpoints) /* + KG * K*/;
}

//tex:$$E\left[\Gamma(\vec{x})\right]=\frac{1}{2}k_B\int(C-C_0)^2dS+\frac{1}{2}k_B\int H(C-C_t)\left[(C-C_0)^2-(C_t-C_0)^2\right]dS$$
double capsid_energy_function(const capsid::Harmonics& h)
{
    //Todo: Set Kb + verify implementation
    constexpr double Kb = 1.0;
    const auto size     = h.quadpoints * h.quadpoints;
    double e1_sum = 0.0, e2_sum = 0.0;
    for (u_t_ i = 0; i < size; ++i)
    {
        const auto C    = h.meancurve[i];
        const auto dS   = h.surface_differentials[i];
        const auto cc0  = C - h.C0;
        const auto cct  = C - h.Ct;
        const auto ctc0 = h.Ct - h.C0;
        e1_sum += cc0 * cc0 * dS;
        e2_sum += C * cct * (cc0 * cc0 - ctc0 * ctc0) * dS;
    }
    return capsid::NormalizeIntegral(e1_sum + e2_sum, h.quadpoints) * .5 * Kb;
}
