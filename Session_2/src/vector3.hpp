// Minimal 3D vector type used throughout the Boris pusher: position, velocity,
// and field values are all Vector3, and the pusher itself is expressed purely
// in terms of the operators below (+, -, scalar *, dot, cross, norm).
#pragma once

#include <cmath>

struct Vector3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

inline Vector3 operator+(const Vector3 &a, const Vector3 &b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

inline Vector3 operator-(const Vector3 &a, const Vector3 &b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

inline Vector3 operator*(double scalar, const Vector3 &v) {
    return {scalar * v.x, scalar * v.y, scalar * v.z};
}

inline double dot(const Vector3 &a, const Vector3 &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vector3 cross(const Vector3 &a, const Vector3 &b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline double norm(const Vector3 &v) {
    return std::sqrt(dot(v, v));
}
