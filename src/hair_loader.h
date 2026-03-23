#pragma once
#include "math.h"
#include <string>
#include <vector>

namespace m3hair {

struct HairData {
    std::vector<float4> vertices;   // x, y, z, radius — Embree native layout
    std::vector<int>    indices;    // first CP index of each B-spline segment
    int                 num_curves   = 0;
    int                 num_segments = 0;
};

// Load a .m3hair file. offset is added to every vertex position.
HairData load_m3hair(const std::string& path, vec3 offset = {});

} // namespace m3hair
