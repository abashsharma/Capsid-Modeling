#include <fstream>
#include <iostream>
#include <stdexcept>
#include <functional>

#include "Capsid.h"
#include "Files.h"

namespace
{

void _get_xyz(const capsid::Harmonics& h, const std::function<void(double, double, double)>& writer)
{
    const auto step = Q_STEP(h.quadpoints);

    for (u_t_ k = 0; k < h.quadpoints; ++k)
    {
        const auto theta = capsid::theta_ang(k, step);
        const auto cos_t = std::cos(theta);
        const auto sin_t = std::sin(theta);
        for (u_t_ l = 0; l < h.quadpoints; ++l)
        {
            const auto phi   = capsid::phi_ang(l, step);
            const auto cos_p = std::cos(phi);
            const auto sin_p = std::sin(phi);
            auto r = h.radius[capsid::Index2D(k, l, h.quadpoints)];
            // Spherical to cartesian coordinates.
            //tex:$$ x = r sin(\theta) cos(\phi) \\
            // y = r sin(\phi) sin(\theta) \\
            // z = rcos(\theta)$$
            auto x = r * sin_t * cos_p;
            auto y = r * sin_p * sin_t;
            auto z = r * cos_t;
            writer(x, y, z);
        }
    }

}

void _save_radii_xyz(const capsid::Harmonics& h, const fs::path& ofile, bool app)
{
    std::ofstream o(ofile, (app ? std::ios::app : std::ios::out));
    if (!o.is_open())
    {
        std::cout << "Failed to open file " << ofile << '\n';
        return;
    }
    o << h.radius.size() << "\n\n";
    _get_xyz(h, [&](double x, double y, double z) {
        o << 0 << ' ' << x << ' ' << y << ' ' << z << '\n';
    });
}

void _save_radii_obj(const capsid::Harmonics& h, const fs::path& ofile)
{
    std::ofstream o(ofile);
    if (!o.is_open())
    {
        std::cout << "Failed to open file " << ofile << '\n';
        return;
    }
    _get_xyz(h, [&](double x, double y, double z) {
        o << "v " << x << ' ' << y << ' ' << z << ' ' << 1.0 << '\n';
    });
    // Todo: Figure out what I need to do to ACTUALLY write a useful obj file.
}

}

namespace capsid
{

void SaveRadii(const Harmonics& h, const fs::path& ofile, bool append)
{
    const auto ext = ofile.extension();
    if (ext == ".xyz")
    {
        _save_radii_xyz(h, ofile, append);
    }
    else if (ext == ".obj")
    {
        _save_radii_obj(h, ofile);
    }
    else
    {
        std::cout << "Unrecognized file format\n";
    }
}

const char* SupportedFormats() noexcept
{
    return ".xyz";
}

}
