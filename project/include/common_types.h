#pragma once
#include <vector>
#include <cmath>
#include <limits> // Required for numeric_limits

// Basic type definitions
using Real = double;
constexpr Real REAL_MAX = std::numeric_limits<Real>::max();
constexpr Real REAL_MIN = std::numeric_limits<Real>::lowest(); // Use lowest for negative infinity representation if needed
constexpr Real REAL_EPSILON = std::numeric_limits<Real>::epsilon();

// Simple 2D Vector struct
struct Vector2D {
    Real x = 0.0;
    Real y = 0.0;

    Vector2D& operator+=(const Vector2D& rhs) { x += rhs.x; y += rhs.y; return *this; }
    Vector2D& operator-=(const Vector2D& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
    Vector2D& operator*=(Real scalar) { x *= scalar; y *= scalar; return *this; }
    Vector2D& operator/=(Real scalar) { x /= scalar; y /= scalar; return *this; }
};

inline Vector2D operator+(Vector2D lhs, const Vector2D& rhs) { lhs += rhs; return lhs; }
inline Vector2D operator-(Vector2D lhs, const Vector2D& rhs) { lhs -= rhs; return lhs; }
inline Vector2D operator*(Vector2D lhs, Real scalar) { lhs *= scalar; return lhs; }
inline Vector2D operator*(Real scalar, Vector2D rhs) { rhs *= scalar; return rhs; }
inline Vector2D operator/(Vector2D lhs, Real scalar) { lhs /= scalar; return lhs; }
inline Real dot(const Vector2D& a, const Vector2D& b) { return a.x * b.x + a.y * b.y; }
inline Real magnitude_sq(const Vector2D& v) { return dot(v, v); }
inline Real magnitude(const Vector2D& v) { return std::sqrt(magnitude_sq(v)); }