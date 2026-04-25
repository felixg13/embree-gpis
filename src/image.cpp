#include "image.h"

#include <cmath>
#include <fstream>
#include <spdlog/spdlog.h>
#include <stdexcept>
namespace m3hair {

static float reinhard(float x) {
    return x / (1.F + x);
}
static float gamma22(float x) {
    return std::pow(x, 1.F / 2.2F);
}

void write_ppm(const std::string &path, const Image &img) {
    spdlog::info("Writing PPM: {}", path);
    std::ofstream f(path, std::ios::binary);
    if (!f) {
        throw std::runtime_error("Cannot write: " + path);
}

    f << "P6\n" << img.width << ' ' << img.height << "\n255\n";

    for (const Vec3 &px : img.fb) {
        auto encode = [](float v) -> unsigned char {
            float t = gamma22(reinhard(v));
            t       = t < 0.F ? 0.F : (t > 1.F ? 1.F : t);
            return static_cast<unsigned char>((t * 255.F) + 0.5F);
        };
        unsigned char rgb[3] = {encode(px.x), encode(px.y), encode(px.z)};
        f.write(reinterpret_cast<char *>(rgb), 3);
    }
}

} // namespace m3hair
