#pragma once
#include "bsdf_deon.h"
#include "hair_loader.h"
#include "light.h"
#include "math.h"

#include <cstdint>
#include <embree4/rtcore.h>
#include <vector>

namespace m3hair {

struct RNG {
    uint32_t state;
    explicit RNG(uint32_t seed) : state(seed | 1U) {}

    uint32_t next_u32() {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return state;
    }
    float next_f() { return (next_u32() >> 8) * (1.F / (1U << 24)); }
};

Vec3 bspline_tangent(const HairData &hair, unsigned primID, float u);

Vec3 trace_path(const Ray &ray,
                RTCScene scene,
                const DirectionalLight &light,
                const std::vector<const HairData *> &hairs,
                const std::vector<DeonParams> &params,
                int max_depth,
                RNG &rng);

} // namespace m3hair
