// Standalone BSDF unit tests — compiled as a separate executable.
#include "bsdf_deon.h"
#include <cmath>
#include <cstdio>
#include <spdlog/spdlog.h>

using namespace m3hair;

static constexpr float PI = 3.14159265358979f;

// Build a direction from longitudinal angle theta and azimuthal phi
// relative to tangent T = (0,0,1), with reference perpendicular (1,0,0).
static vec3 dir_from_angles(float theta, float phi) {
    float ct = cosf(theta), st = sinf(theta);
    return { ct*cosf(phi), ct*sinf(phi), st };
}

// -----------------------------------------------------------------------
// Test 1: White furnace — sigma_a=0, integrate f*cos(theta_o) over sphere.
// Expected: result <= 1.0 + eps for all wi.
// -----------------------------------------------------------------------
static bool test_white_furnace() {
    DeonParams p;
    p.sigma_a = {0.f, 0.f, 0.f};   // no absorption
    const vec3 T = {0.f, 0.f, 1.f};

    bool passed = true;
    // Test two incident directions
    for (float theta_i : {0.f, 0.3f, -0.5f}) {
        for (float phi_i : {0.f, 1.0f}) {
            vec3 wi = dir_from_angles(theta_i, phi_i);

            // Stratified grid on the sphere parameterised by (theta_o, phi_o)
            const int NT = 60, NP = 120;
            float integral[3] = {0,0,0};
            float dtheta = PI / NT;
            float dphi   = 2.f*PI / NP;

            for (int it = 0; it < NT; ++it) {
                // theta_o in (-pi/2, pi/2)
                float theta_o = -PI*0.5f + (it + 0.5f)*dtheta;
                float cos_to  = cosf(theta_o);
                for (int ip = 0; ip < NP; ++ip) {
                    float phi_o = (ip + 0.5f)*dphi;
                    vec3 wo = dir_from_angles(theta_o, phi_o);
                    vec3 f  = eval_deon(wi, wo, T, p);
                    // dwo = cos(theta_o) dtheta dphi, multiply by |cos(theta_o)| for rendering eq
                    float weight = cos_to * cos_to * dtheta * dphi;
                    integral[0] += f.x * weight;
                    integral[1] += f.y * weight;
                    integral[2] += f.z * weight;
                }
            }
            for (int c = 0; c < 3; ++c) {
                if (integral[c] > 1.02f) {
                    spdlog::error("  FAIL white_furnace ch={} theta_i={:.2f} phi_i={:.2f}: integral={:.4f}",
                               c, theta_i, phi_i, integral[c]);
                    passed = false;
                }
            }
        }
    }
    return passed;
}

// -----------------------------------------------------------------------
// Test 2: Reciprocity — f(wi,wo) == f(wo,wi) within 1e-3.
// -----------------------------------------------------------------------
static bool test_reciprocity() {
    DeonParams p;
    const vec3 T = {0.f, 0.f, 1.f};
    bool passed = true;

    // Structured sample of pairs
    float thetas[] = { 0.0f, 0.3f, -0.4f, 0.7f };
    float phis[]   = { 0.0f, 1.0f,  2.5f, 4.0f };
    for (float ti : thetas) for (float pi_ : phis)
    for (float to : thetas) for (float po : phis) {
        vec3 wi = dir_from_angles(ti, pi_);
        vec3 wo = dir_from_angles(to, po);
        vec3 f_ab = eval_deon(wi, wo, T, p);
        vec3 f_ba = eval_deon(wo, wi, T, p);
        float err = std::max({fabsf(f_ab.x-f_ba.x),
                              fabsf(f_ab.y-f_ba.y),
                              fabsf(f_ab.z-f_ba.z)});
        if (err > 1e-3f) {
            spdlog::error("  FAIL reciprocity err={:.5f}  ti={:.2f} pi={:.2f} to={:.2f} po={:.2f}",
                       err, ti, pi_, to, po);
            passed = false;
        }
    }
    return passed;
}

// -----------------------------------------------------------------------
// Test 3: R-only lobe peaks near specular cone.
// When beta_TT and beta_TRT are set to near-zero, the R lobe should give
// the highest value for wi and wo close to the specular direction.
// -----------------------------------------------------------------------
static bool test_R_lobe_specular() {
    DeonParams p;
    p.sigma_a = {0.f,0.f,0.f};
    p.beta_TT  = 1e-6f;
    p.beta_TRT = 1e-6f;
    const vec3 T = {0.f, 0.f, 1.f};

    // For pure R, with tilt alpha=0 and beta_n small, the peak should be
    // near phi=0 (same azimuth, reflected longitude).
    vec3 wi = dir_from_angles(0.3f, 0.f);
    // Specular direction: theta_o ~ theta_i, phi_o ~ 0 (R reflects symmetrically)
    vec3 wo_spec  = dir_from_angles(-0.3f, 0.f);
    vec3 wo_off   = dir_from_angles( 0.3f, 2.f);

    vec3 f_spec = eval_deon(wi, wo_spec, T, p);
    vec3 f_off  = eval_deon(wi, wo_off,  T, p);

    bool ok = (f_spec.x > f_off.x) && (f_spec.y > f_off.y) && (f_spec.z > f_off.z);
    if (!ok)
        spdlog::error("  FAIL R_lobe: f_spec=({:.4f},{:.4f},{:.4f}) f_off=({:.4f},{:.4f},{:.4f})",
                   f_spec.x, f_spec.y, f_spec.z, f_off.x, f_off.y, f_off.z);
    return ok;
}

// -----------------------------------------------------------------------
int main() {
    spdlog::set_pattern("[%T.%e] [%^%l%$] %v");
    spdlog::info("=== d'Eon BSDF unit tests ===");

    bool ok1 = test_white_furnace();
    spdlog::log(ok1 ? spdlog::level::info : spdlog::level::err, "[{}] white_furnace",    ok1 ? "PASS" : "FAIL");

    bool ok2 = test_reciprocity();
    spdlog::log(ok2 ? spdlog::level::info : spdlog::level::err, "[{}] reciprocity",      ok2 ? "PASS" : "FAIL");

    bool ok3 = test_R_lobe_specular();
    spdlog::log(ok3 ? spdlog::level::info : spdlog::level::err, "[{}] R_lobe_specular",  ok3 ? "PASS" : "FAIL");

    bool all = ok1 && ok2 && ok3;
    spdlog::log(all ? spdlog::level::info : spdlog::level::err,
                all ? "All tests PASSED." : "Some tests FAILED.");
    return all ? 0 : 1;
}
