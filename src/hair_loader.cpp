#include "hair_loader.h"

#include <cstddef>
#include <cstdio>
#include <fstream>
#include <spdlog/spdlog.h>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace m3hair {

static void flush_block(HairData &out,
                        std::vector<Float4> &block,
                        int &variant,
                        Vec3 offset,
                        float cos_ry,
                        float sin_ry) {
    if (block.empty()) {
        return;
}

    if (variant == 0) {
        variant = static_cast<int>(block.size()) == 4 ? 4 : 67;
        if (variant == 67) {
            out.vertices.reserve(static_cast<std::size_t>(10'000 * 67));
            out.indices.reserve(static_cast<std::size_t>(10'000 * 64));
        } else {
            out.vertices.reserve(static_cast<std::size_t>(1'100'000 * 4));
            out.indices.reserve(1'100'000);
        }
    }

    const int base = static_cast<int>(out.vertices.size());

    for (auto &v : block) {
        float rx = (v.x * cos_ry) + (v.z * sin_ry);
        float rz = (-v.x * sin_ry) + (v.z * cos_ry);
        v.x      = rx + offset.x;
        v.y      = v.y + offset.y;
        v.z      = rz + offset.z;
        out.vertices.push_back(v);
    }

    if (variant == 4) {
        out.indices.push_back(base);
        out.num_segments += 1;
    } else {
        const int n_segs = static_cast<int>(block.size()) - 3;
        for (int i = 0; i < n_segs; ++i) {
            out.indices.push_back(base + i);
}
        out.num_segments += n_segs;
    }

    out.num_curves += 1;
    block.clear();
}

HairData load_m3hair(const std::string &path, Vec3 offset, float rotate_y_deg) {
    std::ifstream f(path);
    if (!f) {
        throw std::runtime_error("Cannot open: " + path);
}

    const float rad    = rotate_y_deg * (std::numbers::pi_v<float> / 180.F);
    const float cos_ry = std::cos(rad);
    const float sin_ry = std::sin(rad);

    HairData out;
    std::vector<Float4> block;
    block.reserve(67);
    int variant = 0;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos) {
            flush_block(out, block, variant, offset, cos_ry, sin_ry);
            continue;
        }
        Float4 v{};
        if (std::sscanf(line.c_str(), "%f %f %f %f", &v.x, &v.y, &v.z, &v.w) == 4) {
            block.push_back(v);
}
    }
    flush_block(out, block, variant, offset, cos_ry, sin_ry);

    spdlog::info(
        "  Loaded '{}': {} curves, {} segments, {} vertices ({:.1f} MB verts + {:.1f} MB idx)",
        path,
        out.num_curves,
        out.num_segments,
        out.vertices.size(),
        out.vertices.size() * sizeof(Float4) / 1e6,
        out.indices.size() * sizeof(int) / 1e6);

    return out;
}

} // namespace m3hair
