#pragma once
#include <cmath>

// ---------------------------------------------------------------------------
// vec2
// ---------------------------------------------------------------------------
struct vec2 {
    float x{}, y{};
    vec2() = default;
    vec2(float x, float y) : x(x), y(y) {}
    vec2 operator+(const vec2 &o) const { return {x + o.x, y + o.y}; }
    vec2 operator-(const vec2 &o) const { return {x - o.x, y - o.y}; }
    vec2 operator*(float t) const { return {x * t, y * t}; }
    vec2 operator-() const { return {-x, -y}; }
};
inline vec2 operator*(float t, const vec2 &v) {
    return v * t;
}
inline float dot(const vec2 &a, const vec2 &b) {
    return a.x * b.x + a.y * b.y;
}
inline float length(const vec2 &v) {
    return std::sqrt(dot(v, v));
}

// ---------------------------------------------------------------------------
// float4 — Embree vertex layout: x, y, z, radius
// ---------------------------------------------------------------------------
struct float4 {
    float x, y, z, w;
};

// ---------------------------------------------------------------------------
// vec3
// ---------------------------------------------------------------------------
struct vec3 {
    float x{}, y{}, z{};

    vec3() = default;
    vec3(float x, float y, float z) : x(x), y(y), z(z) {}
    explicit vec3(float v) : x(v), y(v), z(v) {}

    vec3 operator+(const vec3 &o) const { return {x + o.x, y + o.y, z + o.z}; }
    vec3 operator-(const vec3 &o) const { return {x - o.x, y - o.y, z - o.z}; }
    vec3 operator*(float t) const { return {x * t, y * t, z * t}; }
    vec3 operator*(const vec3 &o) const { return {x * o.x, y * o.y, z * o.z}; }
    vec3 operator/(float t) const { return {x / t, y / t, z / t}; }
    vec3 &operator+=(const vec3 &o) {
        x += o.x;
        y += o.y;
        z += o.z;
        return *this;
    }
    vec3 &operator*=(float t) {
        x *= t;
        y *= t;
        z *= t;
        return *this;
    }
    vec3 operator-() const { return {-x, -y, -z}; }
};

inline vec3 operator*(float t, const vec3 &v) {
    return v * t;
}

inline float dot(const vec3 &a, const vec3 &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline float length(const vec3 &v) {
    return std::sqrt(dot(v, v));
}
inline vec3 normalize(const vec3 &v) {
    return v / length(v);
}
inline vec3 cross(const vec3 &a, const vec3 &b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline vec3 abs(const vec3 &v) {
    return {std::fabs(v.x), std::fabs(v.y), std::fabs(v.z)};
}
inline vec3 clamp(const vec3 &v, float lo, float hi) {
    auto c = [&](float x) {
        return x < lo ? lo : (x > hi ? hi : x);
    };
    return {c(v.x), c(v.y), c(v.z)};
}

// ---------------------------------------------------------------------------
// Ray
// ---------------------------------------------------------------------------
struct Ray {
    vec3 origin;
    vec3 dir;
    float tnear = 1e-4f;
    float tfar  = 1e30f;
};
