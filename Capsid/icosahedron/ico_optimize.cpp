#include "ico_optimize.h"
#include <fstream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>
#include <unordered_set>
#include <stdexcept>

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

//Convert harmonics to xyz
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
void ico_optimize(capsid::Harmonics& h)
{
    
    //Opt Parameters
    constexpr auto NSAMPLES = 1000;
    const auto kT = 0.0005; 
    const auto fname{ "minimize_ico.xyz" };

    auto Kcalc = Calculate_MeanCurve(h);        //Initialize h
    capsid::SaveRadii(h, fname);
    
    std::vector<Vec3> h_xyz = get_xyz(h);        //Get xyz from h
    

    //Generate icosaherdon
    int ico_points = 400;
    auto ico_vertices = generate_icosahedron_surface_points(ico_points);
    //Print out to be sure
    std::ofstream ofs("icoVertices.xyz");
    ofs << ico_vertices.size() << "\nIcoVertices\n";
    for (const auto& p : ico_vertices) {
        ofs << 0 << " " << p.x << " " << p.y << " " << p.z << "\n";
    }
    ofs.close();
    
    // Map all ico_points to nearest h_xyz 
    std::vector<size_t> matched_indices;               // stores matched h_xyz index per ico_vertex
    std::unordered_set<size_t> used_ico_indices;         // track used ico_xyz indices

    for (const auto& v : h_xyz) {
        double min_dist = std::numeric_limits<double>::max();
        size_t best_j = ico_vertices.size();  // invalid default, means no match found

        for (size_t j = 0; j < ico_vertices.size(); ++j) {
            if (used_ico_indices.count(j)) continue;    // skip already matched h_xyz

            double dist = (v - ico_vertices[j]).norm_squared();
            if (dist < min_dist) {
                min_dist = dist;
                best_j = j;
            }
        }

        if (best_j == ico_vertices.size()) {
            // No available unmatched h_xyz points left — skip this ico_vertex
            continue;
        }

        matched_indices.push_back(best_j);
        used_ico_indices.insert(best_j);
    }

    //Just for safety check
    if (matched_indices.empty()) {
        throw std::runtime_error("No ico_vertices matched to h_xyz points.");
    }

    // Reference points per mapping for the ico_vertices so that there is one-one
    std::vector<Vec3> ref_points;
    for (auto idx : matched_indices)
        ref_points.push_back(ico_vertices[idx]);


    //To store previous state
    Vec a0{ std::move(h.a)};
    auto best_rms = compute_rmsd(ref_points, h_xyz);

    std::uniform_real_distribution<double> dAdist(-.01, .01);
    std::uniform_int_distribution<> accept(0, 100);
    std::uniform_int_distribution<> toPerturb(1, h.a.size());

    std::ofstream ofs1("rms.xyz"); //To putput rms later
    ofs1 << std::scientific << std::setprecision(8);
    int accepted = 0;
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

        //Calculate the harmonics again with the new 'a' vector
        Kcalc  = Calculate_MeanCurve(h);
        std::vector<Vec3> new_h_xyz = get_xyz(h);

        auto new_rms = compute_rmsd(ref_points, new_h_xyz);
        
        const auto d_rms = new_rms - best_rms;
        
        // match kT to tolerance
        // tune kt based on typical energy variance as make proposal
        // set kt within 1 stdev of a typical dE
        //Print out to be sure   
        
        if (new_rms < best_rms || accept(capsid::Generator()) / 100.0 < std::exp(-d_rms / kT)) {
            best_rms = new_rms;
            ++accepted;
            ofs1 << accepted<<" "<<best_rms<<"\n";
            capsid::SaveRadii(h, fname, true);
            a0 = std::move(h.a);
        }
    }
    ofs1.close();
    std::cout << "Accepted steps: " << accepted << "/" << NSAMPLES << "\n";
}
