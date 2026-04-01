#pragma once
#include "math.h"

#include <string>
#include <vector>

namespace m3hair {

struct Image {
    int width{}, height{};
    std::vector<vec3> fb;

    Image() = default;
    Image(int w, int h) : width(w), height(h), fb(w * h, vec3(0.f)) {}

    vec3 &at(int x, int y) { return fb[y * width + x]; }
    const vec3 &at(int x, int y) const { return fb[y * width + x]; }
};

void write_ppm(const std::string &path, const Image &img);

} // namespace m3hair
