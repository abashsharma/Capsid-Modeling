#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <iomanip>
#include "harmonics.hpp"
#include "ico_generator.h"
#include "capsid.h"

double compute_rmsd(const std::vector<Vec3D>& ref, const std::vector<Vec3D>& target) {
    double rms = 0.0;
    for (size_t i = 0; i < ref.size(); ++i) {
        Vec3D diff = ref[i] - target[i];
        rms += diff.norm_squared();
    }
    return std::sqrt(rms / ref.size());
}

std::vector<Vec3D> center_points(const std::vector<Vec3D>& points) {
    Vec3D centroid(0, 0, 0);
    for (const auto& p : points)
        centroid = centroid + p;
    centroid = centroid / points.size();

    std::vector<Vec3D> centered;
    for (const auto& p : points)
        centered.push_back(p - centroid);
    return centered;
}

void ico_optimize(capsid::Harmonics& h,
                  const std::vector<Vec3D>& ico_vertices,
                  size_t max_iterations = 100) {
    std::vector<Vec3D> h_xyz = center_points(h.to_xyz());
    std::vector<Vec3D> centered_ico = center_points(ico_vertices);

    // Create fixed nearest-point mapping: ico_vertex[i] → h_xyz[nearest]
    std::vector<size_t> matched_indices;
    for (const auto& v : centered_ico) {
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

    // Prepare reference points from initial mapping
    std::vector<Vec3D> ref_points;
    for (auto idx : matched_indices)
        ref_points.push_back(h_xyz[idx]);

    double best_rms = compute_rmsd(centered_ico, ref_points);

    for (size_t iter = 0; iter < max_iterations; ++iter) {
        // Perturb harmonics (or rely on external update if called inside loop)
        h.update();  // Assumes this modifies h.radius or other internal state

        std::vector<Vec3D> new_h_xyz = center_points(h.to_xyz());

        // Get target points based on original mapping
        std::vector<Vec3D> new_ref;
        for (auto idx : matched_indices)
            new_ref.push_back(new_h_xyz[idx]);

        double new_rms = compute_rmsd(centered_ico, new_ref);

        if (new_rms < best_rms) {
            best_rms = new_rms;

            // Save accepted structure
            std::ofstream ofs("ico_opt_step_" + std::to_string(iter) + ".xyz");
            ofs << new_h_xyz.size() << "\nOptimized step " << iter << "\n";
            for (const auto& p : new_h_xyz)
                ofs << "0 " << std::fixed << std::setprecision(6)
                    << p.x << " " << p.y << " " << p.z << "\n";

            std::cout << "Iteration " << iter << ": RMSD improved to " << best_rms << "\n";
        }
    }
}
