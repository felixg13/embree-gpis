#pragma once
#include "math.h"

namespace m3hair {

struct DirectionalLight {
    vec3 direction; // unit vector TOWARD the light source
    vec3 radiance;  // RGB

    static DirectionalLight make_default() {
        return {normalize({0.f, 1.f, 2.f}), {80.f, 80.f, 80.f}};
    }
};

} // namespace m3hair
