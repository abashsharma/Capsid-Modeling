#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <numbers>
#include <thread>
#include <random>
#include <chrono>
#include <mutex>
#include <fstream>

#include "icosahedron/ico_generator.h"	//generates the icosahedron
#include "icosahedron/ico_optimize.h"	//optimize with ico
#include "Sequence.h"
#include "Utils.h"
#include "Capsid.h"
#include "Files.h"
#include "Optimize.h"



namespace ch = std::chrono;

void usage(const char* invalid);

i_t_ zero_check(const char* name, const char* val);

int main(int argv, char** argc)
{
    double MaxC0 = 10.0;
    double MaxCt = 10.0;
    double CStep = 0.1;
    u_t_ N       = 3; // Default
    u_t_ qpoints = 20;
    std::string method{ "NCG" };

    //Todo: Remove RunMC in the future
    bool RunMC = false;
    bool RunIco=false;

    //Todo: Get a better commandline parser
    for (int arg = 1; arg < argv; ++arg)
    {
        std::string_view opt(argc[arg]);
        if (opt.compare("-N") == 0)
        {
            arg++;
            N = static_cast<u_t_>(zero_check("N", argc[arg]));
        }
        else if (opt.compare("-C0") == 0)
        {
            arg++;
            MaxC0 = std::abs(std::stod(argc[arg]));
        }
        else if (opt.compare("-Ct") == 0)
        {
            arg++;
            MaxCt = std::abs(std::stod(argc[arg]));
        }
        else if (opt.compare("S") == 0)
        {
            arg++;
            CStep = std::abs(std::stod(argc[arg]));
        }
        else if (opt.compare("-Q") == 0)
        {
            arg++;
            qpoints = static_cast<u_t_>(zero_check("Q", argc[arg]));
        }
        else if (opt.compare("-MC") == 0)
        {
            RunMC = true;
        }
        else if (opt.compare("-M") == 0)
        {
            arg++;
            method = argc[arg];
        }
	else if (opt.compare("-ico") == 0)
        {
            RunIco=true;
        }
        else
        {
            usage(argc[arg]);
        }
    }
    // Total number of m and n values in the array
    u_t_ NMode   = capsid::EnergyModes(N);
    u_t_ NPoints = qpoints * qpoints;

    std::cout << "Initial parameters:"
        << "\n    Harmonics degree: " << N
        << "\n    Energy modes: " << NMode
        << "\n    Number of points: " << NPoints
        << '\n';

	
    /*
    //define Capsid object with the give parameters. each a is generated uniformly between 0 and 0.5
    capsid::Harmonics h(qpoints, NMode);
    std::uniform_real_distribution<double> dist(0.0, 0.5);
    for (auto& a_ : h.a)
    {
      a_ = dist(capsid::Generator());
    }
    h.C0 = 1.0;
    h.a[0] = 2.0 * std::sqrt(std::numbers::pi);
    h.a[1] = 2.0;
    */

    //define Capsid object with the give parameters. each a is generated exponentially
    capsid::Harmonics h(qpoints, NMode);
    //auto a_=capsid::randexp(NMode);
    h.a[0] = 2.0 * std::sqrt(std::numbers::pi);
    h.a = capsid::randexp(NMode, h.a[0]);
    h.C0 = MaxC0;
    h.Ct = MaxCt;
    for (int i=0; i< size(h.a); i++)
    {
      std::cout << h.a[i]<<std::endl;
    }

    const auto K = capsid::Calculate_MeanCurve(h);
    capsid::SaveRadii(h, "Initial.xyz");
    
    if (RunMC)
    {
        const fs::path fname{ "minimization_MC.xyz" };
        const auto start{ ch::high_resolution_clock::now() };

        // Randomly generate a coefficients for minimization
        // drawn from a uniform real distribution
        //std::uniform_real_distribution<double> dist(0.0, .1);
        //for (auto& a_ : h.a)
        //{
        //    a_ = dist(capsid::Generator());
        //}
        //h.C0 = 1.0;
        //h.a[0] = 2.0 * std::sqrt(std::numbers::pi);
	
        //h.a[0]=1.0;
        optimize(h, "MC");

        const auto end{ ch::high_resolution_clock::now() };
        std::cout << "Monte Carlo minimization complete. Took: " << ch::duration_cast<ch::milliseconds>(end-start) << '\n';
        std::cout << "Saving results to " << (fs::current_path() / fname).string() << '\n';

        capsid::SaveRadii(h, fname);
    }

    //For Ico
    if (RunIco)
    {
        const fs::path fname{ "minimization_Ico.xyz" };
	const auto start{ ch::high_resolution_clock::now() };

        ico_optimize(h);

        const auto end{ ch::high_resolution_clock::now() };
        std::cout << "Monte Carlo with Ico minimization complete. Took: " << ch::duration_cast<ch::milliseconds>(end-start) << '\n';
	       
        capsid::SaveRadii(h, fname);
    }

   
    std::cout << "Results:"
        << "\n    Analytical Gauss-Bonnet: " << 4 * std::numbers::pi
        << "\n    Calculated Gauss-Bonnet: " << K
        << '\n';


    return 0;
}

void usage(const char* invalid)
{
    std::cout << "Invalid option: " << invalid 
        << R"(Usage:
    -N:  [+Integer]: Degree of harmonics. (default: 9)
    -Q:  [+Integer]: Number of points to compute. (default: 20)
    -S:  [+Double]:  Step size between curvature values. (default: 0.1)
    -M:              Minimization method. (default: NCG)
    -Ct: [+Double]:  Maximum forced curvature value. (default: 10.0)
    -C0: [+Double]:  Maximum buckling value. (default: 10.0)
    -ico	     Minimize with MC icosaherdon
)";
    std::exit(-1);
}

i_t_ zero_check(const char* vname, const char* val)
{
    auto l = std::stol(val);
    if (l <= 0)
    {
        std::cout << vname << " must be greater than 0\n";
        usage(vname);
    }
    return l;
}
