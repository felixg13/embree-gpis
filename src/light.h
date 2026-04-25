#pragma once
#include "math.h"

namespace m3hair {

struct DirectionalLight {
    Vec3 direction;
    Vec3 radiance;

    static DirectionalLight make_default() {
        return {.direction=normalize({0.F, 1.F, 2.F}), .radiance={80.F, 80.F, 80.F}};
    }
};

} // namespace m3hair
