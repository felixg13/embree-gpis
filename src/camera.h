#pragma once
#include "math.h"

#include <cmath>

namespace m3hair {

struct Camera {
    vec3 origin;
    vec3 lower_left;
    vec3 horizontal;
    vec3 vertical;
    int width{}, height{};

    Ray generate_ray(float px, float py) const {
        float u  = px / static_cast<float>(width);
        float v  = (static_cast<float>(height) - 1.0f - py) / static_cast<float>(height);
        vec3 dir = lower_left + u * horizontal + v * vertical - origin;
        return Ray{origin, normalize(dir)};
    }
};

inline Camera make_camera(
    int w, int h, vec3 origin, vec3 target, float fov_deg = 45.f, vec3 up = {0.f, 1.f, 0.f}) {
    const float aspect      = static_cast<float>(w) / static_cast<float>(h);
    const float half_height = std::tan(fov_deg * 0.5f * (3.14159265f / 180.f));
    const float half_width  = aspect * half_height;

    const vec3 forward = normalize(target - origin);
    const vec3 right   = normalize(cross(forward, up));
    const vec3 up_corr = cross(right, forward);

    Camera cam;
    cam.origin     = origin;
    cam.width      = w;
    cam.height     = h;
    cam.horizontal = right * (2.f * half_width);
    cam.vertical   = up_corr * (2.f * half_height);
    cam.lower_left = origin + forward - right * half_width - up_corr * half_height;
    return cam;
}

inline Camera make_default_camera(int w, int h) {
    return make_camera(w, h, {0.f, 5.f, 20.f}, {0.f, 5.f, 0.f}, 45.f);
}

} // namespace m3hair
