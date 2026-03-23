#include <embree4/rtcore.h>
#include <print>
#include <atomic>
#include <vector>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

#include "bsdf_deon.h"
#include "camera.h"
#include "hair_loader.h"
#include "image.h"
#include "integrator.h"
#include "light.h"
#include "math.h"
#include "scene.h"

int main()
{
    const int W         = 1200, H = 600;
    const int SPP       = 16;
    const int MAX_DEPTH = 4;

    RTCDevice device = rtcNewDevice(nullptr);
    if (!device) { std::print(stderr, "Embree device failed\n"); return 1; }

    // placed the assets manually to have a side by side view
    std::print("Loading assets ...\n");
    m3hair::HairData curl = m3hair::load_m3hair("assets/curl.m3hair", {15.f, 3.f, 6.5f});
    m3hair::HairData hair = m3hair::load_m3hair("assets/hair.m3hair", { 0.f, 0.f, 0.f});

    std::print("Building scene ...\n");
    m3hair::Scene scene(device);
    scene.add_hair(curl);   // geomID 0
    scene.add_hair(hair);   // geomID 1
    scene.commit();

    m3hair::Camera cam = m3hair::make_camera(W, H,
        {7.f, 9.f, 12.f}, {7.f, 9.f, 0.f}, 70.f);

    m3hair::Image img(W, H);
    m3hair::DirectionalLight light = m3hair::DirectionalLight::make_default();

    // One DeonParams per geometry (golden for both for now)
    m3hair::DeonParams p_curl, p_hair;
    std::vector<const m3hair::HairData*> hairs = { &curl, &hair };
    std::vector<m3hair::DeonParams>      dparams = { p_curl, p_hair };

    const int total_tiles = ((W+15)/16) * ((H+15)/16);
    std::atomic<int> done_tiles{0};
    std::print("Rendering {}x{} @ {}spp ...\n", W, H, SPP);

    tbb::parallel_for(
        tbb::blocked_range2d<int>(0, H, 16, 0, W, 16),
        [&](const tbb::blocked_range2d<int>& r) {
            for (int y = r.rows().begin(); y < r.rows().end(); ++y) {
                for (int x = r.cols().begin(); x < r.cols().end(); ++x) {
                    m3hair::RNG rng(y * W + x + 1);
                    vec3 acc(0.f);
                    for (int s = 0; s < SPP; ++s) {
                        Ray ray = cam.generate_ray(x + rng.next_f(), y + rng.next_f());
                        acc = acc + m3hair::trace_path(ray, scene.handle(), light,
                                                       hairs, dparams, MAX_DEPTH, rng);
                    }
                    img.at(x, y) = acc * (1.f / SPP);
                }
            }
            int t = ++done_tiles;
            if (t % std::max(1, total_tiles/100) == 0)
                std::print(stderr, "\r  {:.0f}%   ", 100.f*t/total_tiles);
        });
    std::print(stderr, "\r  100%   \n");

    m3hair::write_ppm("out.ppm", img);
    std::print("Wrote out.ppm\n");

    rtcReleaseDevice(device);
    return 0;
}
