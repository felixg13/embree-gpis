#pragma once
#include <cmath>

struct Vec2 {
    float x{}, y{};
    Vec2() = default;
    Vec2(float x, float y) : x(x), y(y) {}
    Vec2 operator+(const Vec2 &o) const { return {x + o.x, y + o.y}; }
    Vec2 operator-(const Vec2 &o) const { return {x - o.x, y - o.y}; }
    Vec2 operator*(float t) const { return {x * t, y * t}; }
    Vec2 operator-() const { return {-x, -y}; }
};
inline Vec2 operator*(float t, const Vec2 &v) {
    return v * t;
}
inline float dot(const Vec2 &a, const Vec2 &b) {
    return (a.x * b.x) + (a.y * b.y);
}
inline float length(const Vec2 &v) {
    return std::sqrt(dot(v, v));
}

// x, y, z, radius — Embree vertex layout
struct Float4 {
    float x, y, z, w;
};

struct Vec3 {
    float x{}, y{}, z{};

    Vec3() = default;
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
    explicit Vec3(float v) : x(v), y(v), z(v) {}

    Vec3 operator+(const Vec3 &o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3 &o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(float t) const { return {x * t, y * t, z * t}; }
    Vec3 operator*(const Vec3 &o) const { return {x * o.x, y * o.y, z * o.z}; }
    Vec3 operator/(float t) const { return {x / t, y / t, z / t}; }
    Vec3 &operator+=(const Vec3 &o) {
        x += o.x;
        y += o.y;
        z += o.z;
        return *this;
    }
    Vec3 &operator*=(float t) {
        x *= t;
        y *= t;
        z *= t;
        return *this;
    }
    Vec3 operator-() const { return {-x, -y, -z}; }
};

inline Vec3 operator*(float t, const Vec3 &v) {
    return v * t;
}
inline float dot(const Vec3 &a, const Vec3 &b) {
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}
inline float length(const Vec3 &v) {
    return std::sqrt(dot(v, v));
}
inline Vec3 normalize(const Vec3 &v) {
    return v / length(v);
}
inline Vec3 cross(const Vec3 &a, const Vec3 &b) {
    return {(a.y * b.z) - (a.z * b.y), (a.z * b.x) - (a.x * b.z), (a.x * b.y) - (a.y * b.x)};
}
inline Vec3 abs(const Vec3 &v) {
    return {std::fabs(v.x), std::fabs(v.y), std::fabs(v.z)};
}
inline Vec3 clamp(const Vec3 &v, float lo, float hi) {
    auto c = [&](float x) {
        if (x < lo) {
            return lo;
        }
        if (x > hi) {
            return hi;
        }
        return x;
    };
    return {c(v.x), c(v.y), c(v.z)};
}

struct Ray {
    Vec3 origin;
    Vec3 dir;
    float tnear = 1e-4F;
    float tfar  = 1e30F;
};
