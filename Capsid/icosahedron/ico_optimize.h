#pragma once

#ifndef ICO_OPTIMIZE_HPP
#define ICO_OPTIMIZE_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <iomanip>
#include "ico_generator.h"
#include "../Capsid.h"
#include "../Files.h"


// Compute RMSD between two vectors of Vec3D points
double compute_rmsd(const std::vector<Vec3>& ref, const std::vector<Vec3>& target);

// Center points by subtracting centroid from all points
std::vector<Vec3> center_points(const std::vector<Vec3>& points);

// Optimize the harmonics shape to match ico vertices via fixed nearest-point assignment
// h: harmonic shape to optimize
// ico_vertices: fixed icosahedron vertex points
// max_iterations: number of optimization iterations (default 100)
void ico_optimize(capsid::Harmonics& h);

#endif // ICO_OPTIMIZE_HPP
