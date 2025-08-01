#include "ico_optimize.h"
#include "Capsid.h"
#include <fstream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>

// Compute RMSD between two point sets
double compute_rmsd(const std::vector<Vec3>& ref, const std::vector<Vec3>& target) {
    double rms = 0.0;
    //For safety, insclude a check
    if (ref.size() != target.size()) {
    throw std::invalid_argument("RMSD vectors must be same length");
    }
    for (size_t i = 0; i < ref.size(); ++i) {
        Vec3 diff = ref[i] - target[i];
        rms += diff.norm_squared();
    }
    return std::sqrt(rms / ref.size());
}

std::vector<Vec3> get_xyz(const capsid::Harmonics& h)
{
    std::vector<Vec3> h_points(h.quadpoints * h.quadpoints);
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

            const auto idx = capsid::Index2D(k, l, h.quadpoints);
            auto r = h.radius[idx];
            h_points[idx].x = r * sin_t * cos_p;
            h_points[idx].y = r * sin_p * sin_t;
            h_points[idx].z = r * cos_t;
        }
    }
    return h_points;
}

/*
// Center points by subtracting centroid
std::vector<Vec3> center_points(const std::vector<Vec3>& points) {
    Vec3 centroid(0, 0, 0);
    for (const auto& p : points)
        centroid = centroid + p;
    centroid = centroid * (1.0 / points.size());

    std::vector<Vec3> centered;
    for (const auto& p : points)
        centered.push_back(p - centroid);
    return centered;
}
*/

// Optimization function
void ico_optimize(capsid::Harmonics& h, const std::vector<Vec3>& ico_vertices)
{
    std::vector<Vec3> h_xyz = get_xyz(h);
    //std::vector<Vec3> centered_ico = center_points(ico_vertices);

    // Map ico vertex to nearest h_xyz point index
    std::vector<size_t> matched_indices;
    for (const auto& v : ico_vertices) {
        double min_dist = std::numeric_limits<double>::max();
        size_t best_j = 0;
        for (size_t j = 0; j < h_xyz.size(); ++j) {
            double dist = (v - h_xyz[j]).norm_squared();
            if (dist < min_dist) {
                min_dist = dist;
                best_j = j;
            }
        }
        matched_indices.push_back(best_j);
    }

    // Reference points per mapping
    std::vector<Vec3> ref_points;
    for (auto idx : matched_indices)
        ref_points.push_back(h_xyz[idx]);


    //dist Dwrite(std::format("d_{}_{}.csv", h.C0, h.Ct));
    const auto fname{ "minimize_ico.xyz" };

        
    constexpr auto NSAMPLES = 1000;
    const auto kT = 0.0005; 
    std::uniform_real_distribution<double> dAdist(-.01, .01);
    std::uniform_int_distribution<> accept(0, 100);
    std::uniform_int_distribution<> toPerturb(1, h.a.size());

    auto rms_calc = compute_rmsd(ico_vertices, ref_points);
    
    capsid::SaveRadii(h, fname);
    Vec a0{ std::move(h.a)};
    double best_rms = compute_rmsd(ico_vertices, ref_points);

    for (int i = 0; i < NSAMPLES; ++i)
    {
        int ntb = toPerturb(capsid::Generator());
        Vec anew{a0};
        for (int j = 0; j < ntb; ++j)
        {
            auto idx   = toPerturb(capsid::Generator()) - 1;
            anew[idx] += dAdist(capsid::Generator());
        }
        // set a subset of as from a distribution
        h.a = std::move(anew);
    

        std::vector<Vec3> new_h_xyz = get_xyz(h);

        std::vector<Vec3> new_ref;
        for (auto idx : matched_indices)
            new_ref.push_back(new_h_xyz[idx]);

        double new_rms = compute_rmsd(ico_vertices, new_ref);

        const auto d_rms = new_new - best_rms;
        // match kT to tolerance
        // tune kt based on typical energy variance as make proposal
        // set kt within 1 stdev of a typical dE
        if (new_rms < best_rms || accept(capsid::Generator()) / 100.0 < std::exp(-d_rms / kT)) {
            best_rms = new_rms;
            capsid::SaveRadii(h, fname, true);
            a0 = std::move(h.a);
        }
    }
}
