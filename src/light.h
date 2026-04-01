#pragma once
#include "math.h"

namespace m3hair {

struct DirectionalLight {
    vec3 direction;
    vec3 radiance;

    static DirectionalLight make_default() {
        return {normalize({0.f, 1.f, 2.f}), {80.f, 80.f, 80.f}};
    }
};

} // namespace m3hair
