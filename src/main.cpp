#include <embree4/rtcore.h>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <print>
#include <string>
#include <string_view>
#include <vector>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>
#include <spdlog/spdlog.h>

#include "bsdf_deon.h"
#include "camera.h"
#include "user_geo.h"
#include "hair_loader.h"
#include "image.h"
#include "integrator.h"
#include "light.h"
#include "math.h"
#include "scene.h"

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------
struct Args {
    std::vector<std::string> files;
    std::vector<float>       rotations;  // per-file Y-axis rotation in degrees
    int   width   = 1200;
    int   height  =  600;
    int   spp     =   16;
    int   depth   =    4;
    float spacing =  2.0f;
    std::string output = "renders/out.ppm";
    std::string mode   = "embree";
};

static void usage(const char* argv0) {
    std::print("Usage: {} [OPTIONS] file1.m3hair [file2.m3hair ...]\n\n"
               "  --width  <int>     image width  (default: 1200)\n"
               "  --height <int>     image height (default: 600)\n"
               "  --spp    <int>     samples per pixel (default: 16)\n"
               "  --depth  <int>     max path depth (default: 4)\n"
               "  --spacing <float>  X spacing between assets (default: 2.0)\n"
               "  --rotate  IDX:DEG  Y-axis rotation for mesh at index IDX (e.g. 1:90)\n"
               "  -o <path>          output PPM (default: renders/out.ppm)\n"
               "  --mode <string>    embree|raymarching (default: embree)\n",
               argv0);
}

static Args parse_args(int argc, char* argv[]) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string_view s = argv[i];
        auto next = [&]() -> std::string_view {
            if (i+1 >= argc) { spdlog::error("Missing value for {}", s); std::exit(1); }
            return argv[++i];
        };
        if      (s == "--width")   a.width   = std::stoi(std::string(next()));
        else if (s == "--height")  a.height  = std::stoi(std::string(next()));
        else if (s == "--spp")     a.spp     = std::stoi(std::string(next()));
        else if (s == "--depth")   a.depth   = std::stoi(std::string(next()));
        else if (s == "--spacing") a.spacing = std::stof(std::string(next()));
        else if (s == "-o")        a.output  = std::string(next());
        else if (s == "--mode")    a.mode    = std::string(next());
        else if (s == "--rotate") {
            std::string val = std::string(next());
            auto colon = val.find(':');
            if (colon == std::string::npos) {
                spdlog::error("--rotate requires IDX:DEG format"); std::exit(1);
            }
            int idx = std::stoi(val.substr(0, colon));
            float deg = std::stof(val.substr(colon + 1));
            if (idx >= (int)a.rotations.size()) a.rotations.resize(idx + 1, 0.f);
            a.rotations[idx] = deg;
        }
        else if (s == "--help" || s == "-h") { usage(argv[0]); std::exit(0); }
        else if (s.starts_with("--")) {
            spdlog::error("Unknown option: {}", s); std::exit(1);
        } else {
            a.files.push_back(std::string(s));
        }
    }
    if (a.files.empty()) {
        // default demo: both assets with manual placement
        a.files = { "assets/curl.m3hair", "assets/hair.m3hair" };
    }
    return a;
}

// ---------------------------------------------------------------------------
// Auto-framing camera — fits all loaded geometry into the frame
// ---------------------------------------------------------------------------
static m3hair::Camera auto_camera(const std::vector<m3hair::HairData>& all_hair,
                                  int W, int H)
{
    float xmin =  1e30f, xmax = -1e30f;
    float ymin =  1e30f, ymax = -1e30f;
    float zmin =  1e30f, zmax = -1e30f;
    for (const auto& h : all_hair) {
        for (const auto& v : h.vertices) {
            xmin = std::min(xmin, v.x); xmax = std::max(xmax, v.x);
            ymin = std::min(ymin, v.y); ymax = std::max(ymax, v.y);
            zmin = std::min(zmin, v.z); zmax = std::max(zmax, v.z);
        }
    }
    float cx = (xmin+xmax)*0.5f, cy = (ymin+ymax)*0.5f;
    float scene_w = xmax - xmin, scene_h = ymax - ymin;
    float aspect  = (float)W / (float)H;
    float fov_deg = 65.f;
    float half_fov = fov_deg * 0.5f * (3.14159265f/180.f);
    // Distance so the wider dimension fills ~85% of the frame
    float d_from_w = (scene_w / aspect) / (2.f * std::tan(half_fov)) / 0.85f;
    float d_from_h = scene_h              / (2.f * std::tan(half_fov)) / 0.85f;
    float d = std::max(d_from_w, d_from_h);
    float cam_z = zmax + d;
    return m3hair::make_camera(W, H,
        {cx, cy, cam_z}, {cx, cy, (zmin+zmax)*0.5f}, fov_deg);
}

// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    spdlog::set_pattern("[%T.%e] [%^%l%$] %v");

    Args a = parse_args(argc, argv);

    RTCDevice device = rtcNewDevice(nullptr);
    if (!device) { spdlog::error("Embree device failed"); return 1; }

    // Load each file with X offset i*spacing
    spdlog::info("Loading {} asset(s)", a.files.size());
    std::vector<m3hair::HairData> all_hair;
    all_hair.reserve(a.files.size());
    for (int i = 0; i < (int)a.files.size(); ++i) {
        vec3  offset  = { i * a.spacing, 0.f, 0.f };
        float rot_deg = i < (int)a.rotations.size() ? a.rotations[i] : 0.f;
        all_hair.push_back(m3hair::load_m3hair(a.files[i], offset, rot_deg));
    }

    spdlog::info("Building scene (mode={})", a.mode);
    m3hair::Scene scene(device);
    for (const auto& h : all_hair) {
        if (a.mode == "raymarching")
            scene.add_user_hair(h);
        else
            scene.add_hair(h);
    }
    scene.commit();
    spdlog::info("Scene committed");

    // Camera: auto-frame all geometry
    m3hair::Camera cam = auto_camera(all_hair, a.width, a.height);

    m3hair::Image img(a.width, a.height);
    m3hair::DirectionalLight light = m3hair::DirectionalLight::make_default();

    // One default DeonParams per geometry
    std::vector<const m3hair::HairData*> hair_ptrs;
    std::vector<m3hair::DeonParams>      dparams;
    for (const auto& h : all_hair) {
        hair_ptrs.push_back(&h);
        dparams.emplace_back();
    }

    const int total_tiles = ((a.width+15)/16) * ((a.height+15)/16);
    std::atomic<int> done_tiles{0};
    spdlog::info("Rendering {}x{} @ {} spp  mode={}  ({} tiles)",
                 a.width, a.height, a.spp, a.mode, total_tiles);

    auto t_start = std::chrono::steady_clock::now();

    tbb::parallel_for(
        tbb::blocked_range2d<int>(0, a.height, 16, 0, a.width, 16),
        [&](const tbb::blocked_range2d<int>& r) {
            for (int y = r.rows().begin(); y < r.rows().end(); ++y) {
                for (int x = r.cols().begin(); x < r.cols().end(); ++x) {
                    m3hair::RNG rng((uint32_t)(y * a.width + x + 1));
                    vec3 acc(0.f);
                    for (int s = 0; s < a.spp; ++s) {
                        Ray ray = cam.generate_ray(x + rng.next_f(), y + rng.next_f());
                        acc = acc + m3hair::trace_path(ray, scene.handle(), light,
                                                       hair_ptrs, dparams,
                                                       a.depth, rng);
                    }
                    img.at(x, y) = acc * (1.f / a.spp);
                }
            }
            int t = ++done_tiles;
            // Log at each 10% milestone — fires exactly once per crossing
            int new_pct = (int)(100.f * t / total_tiles);
            int old_pct = (int)(100.f * (t - 1) / total_tiles);
            if ((new_pct / 10) > (old_pct / 10)) {
                float elapsed = std::chrono::duration<float>(
                    std::chrono::steady_clock::now() - t_start).count();
                spdlog::info("  {:3d}%  ({:.1f}s elapsed)", (new_pct / 10) * 10, elapsed);
            }
        });

    float total_s = std::chrono::duration<float>(
        std::chrono::steady_clock::now() - t_start).count();
    spdlog::info("Render complete in {:.2f}s", total_s);

    std::filesystem::create_directories(
        std::filesystem::path(a.output).parent_path().empty()
            ? std::filesystem::path(".")
            : std::filesystem::path(a.output).parent_path());
    m3hair::write_ppm(a.output, img);
    spdlog::info("Wrote {}", a.output);

    rtcReleaseDevice(device);
    return 0;
}
