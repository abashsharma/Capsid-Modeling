#include "ico_generator.h"
##include <cmath>
##include <vector>
#include <cstdlib>

Vec3::Vec3() : x(0), y(0), z(0) {}
Vec3::Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

Vec3 Vec3::operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
Vec3 Vec3::operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
Vec3 Vec3::operator*(double s) const { return Vec3(x * s, y * s, z * s); }

bool Vec3::operator<(const Vec3& other) const {
    return std::tie(x, y, z) < std::tie(other.x, other.y, other.z);
}

void Vec3::normalize() {
    double len = std::sqrt(x * x + y * y + z * z);
    if (len > 0.0) {
        x /= len;
        y /= len;
        z /= len;
    }
}

static std::vector<Vec3> generate_icosahedron_vertices() {
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    const double a = 1.0;
    const double b = 1.0 / phi;

    std::vector<Vec3> base = {
        { 0,  b, -a}, { b,  a, 0}, {-b,  a, 0}, { 0,  b,  a},
        { 0, -b,  a}, {-a, 0,  b}, { 0, -b, -a}, { a, 0, -b},
        { a, 0,  b}, {-a, 0, -b}, { b, -a, 0}, {-b, -a, 0}
    };

    //  Scale vertices so icosahedron has volume ≈ unit sphere
    const double scale = 1.051462224; // // cube root of (4/3 * pi / (5/12)(3+√5))
    for (auto& v : base) {
        v = v * scale;
    }
  
    return base;   // true icosahedron
}

static std::vector<std::tuple<int, int, int>> generate_icosahedron_faces() {
    return {
        {0, 1, 2}, {3, 2, 1}, {3, 4, 5}, {3, 8, 4}, {0, 6, 7},
        {0, 9, 6}, {4, 10, 11}, {6, 11, 10}, {2, 5, 9}, {11, 9, 5},
        {1, 7, 8}, {10, 8, 7}, {3, 5, 2}, {3, 1, 8}, {0, 2, 9},
        {0, 7, 1}, {6, 9, 11}, {6, 10, 7}, {4, 11, 5}, {4, 8, 10}
    };
}

std::vector<Vec3> generate_icosahedron_surface_points(int desired_points) {
    auto verts = generate_icosahedron_vertices();
    auto faces = generate_icosahedron_faces();
    std::vector<Vec3> points;

    // Estimate subdivision level for ~desired_points total
    int subdivisions = 1;
    while (20 * ((subdivisions + 1) * (subdivisions + 2)) / 2 < desired_points) {
        ++subdivisions;
    }

    for (const auto& [i1, i2, i3] : faces) {
        Vec3 v1 = verts[i1], v2 = verts[i2], v3 = verts[i3];

        for (int i = 0; i <= subdivisions; ++i) {
            for (int j = 0; j <= subdivisions - i; ++j) {
                int k = subdivisions - i - j;

                       double a = double(i) / subdivisions;
                       double b = double(j) / subdivisions;
                       double c = double(k) / subdivisions;
                       Vec3 point = v1 * a + v2 * b + v3 * c;
                       points.push_back(point);
                }
            }
    }

    return points;
}
