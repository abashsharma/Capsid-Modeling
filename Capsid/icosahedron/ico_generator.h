#pragma once

#ifndef ICOSAHEDRON_H
#define ICOSAHEDRON_H

#include <vector>
#include <tuple>

struct Vec3 {
    double x, y, z;

    Vec3();
    Vec3(double x_, double y_, double z_);

    Vec3 operator+(const Vec3& v) const;
    Vec3 operator-(const Vec3& v) const;
    Vec3 operator*(double s) const;

    bool operator<(const Vec3& other) const;

    void normalize();
    double norm_squared() const;
}
};

// Generates approximately `desired_points` distributed across the surface of a geometric icosahedron
std::vector<Vec3> generate_icosahedron_surface_points(int desired_points = 400);

#endif  // ICOSAHEDRON_H
