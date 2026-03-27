#pragma once
#include "math.h"
#include <cstdint>
#include <string>
#include <vector>

namespace m3hair {

struct HairData {
    std::vector<float4> vertices;   // x, y, z, radius — Embree native layout
    std::vector<int>    indices;    // first CP index of each B-spline segment
    int                 num_curves   = 0;
    int                 num_segments = 0;
    float               amplitude    = 0.f;         // noise amplitude (replaces SC_AMPLITUDE)
    uint32_t            seed         = 0xdeadbeef;  // global noise seed (replaces SC_GLOBAL_SEED)
};

// Load a .m3hair file. offset is added to every vertex position.
// rotate_y_deg rotates vertices around the Y-axis before applying offset.
HairData load_m3hair(const std::string& path, vec3 offset = {}, float rotate_y_deg = 0.f);

} // namespace m3hair
