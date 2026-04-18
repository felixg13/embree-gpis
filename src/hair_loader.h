#pragma once
#include "math.h"

#include <cstdint>
#include <string>
#include <vector>

namespace m3hair {

struct HairData {
    std::vector<float4> vertices;
    std::vector<int> indices;
    int num_curves   = 0;
    int num_segments = 0;
    float amplitude  = 0.f;
    float cell_size  = 0.05f;
    uint32_t seed    = 0xdeadbeef;
};

HairData load_m3hair(const std::string &path, vec3 offset = {}, float rotate_y_deg = 0.f);

} // namespace m3hair
