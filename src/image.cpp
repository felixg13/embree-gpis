#include "image.h"
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <spdlog/spdlog.h>
namespace m3hair {

static float reinhard(float x) { return x / (1.f + x); }
static float gamma22(float x) { return std::pow(x, 1.f / 2.2f); }

void write_ppm(const std::string &path, const Image &img) {
  spdlog::info("Writing PPM: {}", path);
  std::ofstream f(path, std::ios::binary);
  if (!f)
    throw std::runtime_error("Cannot write: " + path);

  f << "P6\n" << img.width << ' ' << img.height << "\n255\n";

  for (const vec3 &px : img.fb) {
    auto encode = [](float v) -> unsigned char {
      float t = gamma22(reinhard(v));
      t = t < 0.f ? 0.f : (t > 1.f ? 1.f : t);
      return static_cast<unsigned char>(t * 255.f + 0.5f);
    };
    unsigned char rgb[3] = {encode(px.x), encode(px.y), encode(px.z)};
    f.write(reinterpret_cast<char *>(rgb), 3);
  }
}

} // namespace m3hair
