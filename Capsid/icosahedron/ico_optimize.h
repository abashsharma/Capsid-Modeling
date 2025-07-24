#ifndef ICO_OPTIMIZE_HPP
#define ICO_OPTIMIZE_HPP

#include <vector>
#include "ico_generator.h"    // If needed for Vec3D or related types
#include "capsid.h"           // For capsid related definitions


// Compute RMSD between two vectors of Vec3D points
double compute_rmsd(const std::vector<Vec3D>& ref, const std::vector<Vec3D>& target);

// Center points by subtracting centroid from all points
std::vector<Vec3D> center_points(const std::vector<Vec3D>& points);

// Optimize the harmonics shape to match ico vertices via fixed nearest-point assignment
// h: harmonic shape to optimize
// ico_vertices: fixed icosahedron vertex points
// max_iterations: number of optimization iterations (default 100)
void ico_optimize(capsid::Harmonics& h,
                  const std::vector<Vec3D>& ico_vertices,
                  size_t max_iterations = 100);

#endif // ICO_OPTIMIZE_HPP
