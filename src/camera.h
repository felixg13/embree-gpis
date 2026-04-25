#pragma once
#include "math.h"

#include <cmath>
#include <numbers>

namespace m3hair {

struct Camera {
    Vec3 origin;
    Vec3 lower_left;
    Vec3 horizontal;
    Vec3 vertical;
    int width{}, height{};

    [[nodiscard]] Ray generate_ray(float px, float py) const {
        float u  = px / static_cast<float>(width);
        float v  = (static_cast<float>(height) - 1.0F - py) / static_cast<float>(height);
        Vec3 dir = lower_left + u * horizontal + v * vertical - origin;
        return Ray{.origin=origin, .dir=normalize(dir)};
    }
};

inline Camera make_camera(
    int w, int h, Vec3 origin, Vec3 target, float fov_deg = 45.F, Vec3 up = {0.F, 1.F, 0.F}) {
    const float aspect      = static_cast<float>(w) / static_cast<float>(h);
    const float half_height = std::tan(fov_deg * 0.5F * (std::numbers::pi_v<float> / 180.F));
    const float half_width  = aspect * half_height;

    const Vec3 forward = normalize(target - origin);
    const Vec3 right   = normalize(cross(forward, up));
    const Vec3 up_corr = cross(right, forward);

    Camera cam;
    cam.origin     = origin;
    cam.width      = w;
    cam.height     = h;
    cam.horizontal = right * (2.F * half_width);
    cam.vertical   = up_corr * (2.F * half_height);
    cam.lower_left = origin + forward - right * half_width - up_corr * half_height;
    return cam;
}

inline Camera make_default_camera(int w, int h) {
    return make_camera(w, h, {0.F, 5.F, 20.F}, {0.F, 5.F, 0.F}, 45.F);
}

} // namespace m3hair
