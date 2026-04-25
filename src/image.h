#pragma once
#include "math.h"

#include <cstddef>
#include <string>
#include <vector>

namespace m3hair {

struct Image {
    int width{}, height{};
    std::vector<Vec3> fb;

    Image() = default;
    Image(int w, int h) : width(w), height(h), fb(static_cast<std::size_t>(w * h), Vec3(0.F)) {}

    Vec3 &at(int x, int y) { return fb[(y * width) + x]; }
    [[nodiscard]] const Vec3 &at(int x, int y) const { return fb[(y * width) + x]; }
};

void write_ppm(const std::string &path, const Image &img);

} // namespace m3hair
