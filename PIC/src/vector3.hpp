/*
Pablo Requeijo, September 24 2026.
Session 2, HPC coursework: Boris pusher.
Minimal 3D vector type: position, velocity and field values are all
Vector3, and the pusher is expressed purely in terms of the operators
below. Some (-, norm) are not needed by the pusher but are kept for
convenience.
*/
#pragma once

#include <cmath>

struct Vector3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

// Component-wise sum.
inline Vector3 operator+(const Vector3 &a, const Vector3 &b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

// Component-wise difference.
inline Vector3 operator-(const Vector3 &a, const Vector3 &b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

// Scales every component by a scalar.
inline Vector3 operator*(double scalar, const Vector3 &v) {
    return {scalar * v.x, scalar * v.y, scalar * v.z};
}

// Dot product a . b.
inline double dot(const Vector3 &a, const Vector3 &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Cross product a x b.
inline Vector3 cross(const Vector3 &a, const Vector3 &b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

// Euclidean length |v|.
inline double norm(const Vector3 &v) {
    return std::sqrt(dot(v, v));
}
