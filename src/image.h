#pragma once
#include "math.h"

#include <string>
#include <vector>

namespace m3hair {

struct Image {
    int width{}, height{};
    std::vector<vec3> fb; // linear-light RGB, one vec3 per pixel (row-major)

    Image() = default;
    Image(int w, int h) : width(w), height(h), fb(w * h, vec3(0.f)) {}

    vec3 &at(int x, int y) { return fb[y * width + x]; }
    const vec3 &at(int x, int y) const { return fb[y * width + x]; }
};

// Reinhard tone-map + gamma 2.2, write P6 PPM.
void write_ppm(const std::string &path, const Image &img);

} // namespace m3hair
