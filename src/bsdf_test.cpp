#include "bsdf_deon.h"

#include <cmath>
#include <cstdio>
#include <numbers>
#include <spdlog/spdlog.h>

using namespace m3hair;

static constexpr float PI = std::numbers::pi_v<float>;

static Vec3 dir_from_angles(float theta, float phi) {
    float ct = cosf(theta);
    float st = sinf(theta);
    return {ct * cosf(phi), ct * sinf(phi), st};
}

static bool test_white_furnace() {
    DeonParams p;
    p.sigma_a    = {0.F, 0.F, 0.F};
    const Vec3 t = {0.F, 0.F, 1.F};

    bool passed = true;
    for (float theta_i : {0.F, 0.3F, -0.5F}) {
        for (float phi_i : {0.F, 1.0F}) {
            Vec3 wi = dir_from_angles(theta_i, phi_i);

            const int NT = 60;
            const int NP = 120;
            float integral[3] = {0, 0, 0};
            float dtheta      = PI / NT;
            float dphi        = 2.F * PI / NP;

            for (int it = 0; it < NT; ++it) {
                float theta_o = (-PI * 0.5F) + ((it + 0.5F) * dtheta);
                float cos_to  = cosf(theta_o);
                for (int ip = 0; ip < NP; ++ip) {
                    float phi_o  = (ip + 0.5F) * dphi;
                    Vec3 wo      = dir_from_angles(theta_o, phi_o);
                    Vec3 f       = eval_deon(wi, wo, t, p);
                    float weight = cos_to * cos_to * dtheta * dphi;
                    integral[0] += f.x * weight;
                    integral[1] += f.y * weight;
                    integral[2] += f.z * weight;
                }
            }
            for (int c = 0; c < 3; ++c) {
                if (integral[c] > 1.02F) {
                    spdlog::error(
                        "  FAIL white_furnace ch={} theta_i={:.2f} phi_i={:.2f}: integral={:.4f}",
                        c,
                        theta_i,
                        phi_i,
                        integral[c]);
                    passed = false;
                }
            }
        }
    }
    return passed;
}

static bool test_reciprocity() {
    DeonParams p;
    const Vec3 t = {0.F, 0.F, 1.F};
    bool passed  = true;

    float thetas[] = {0.0F, 0.3F, -0.4F, 0.7F};
    float phis[]   = {0.0F, 1.0F, 2.5F, 4.0F};
    for (float ti : thetas) {
        for (float pi : phis) {
            for (float to : thetas) {
                for (float po : phis) {
                    Vec3 wi   = dir_from_angles(ti, pi);
                    Vec3 wo   = dir_from_angles(to, po);
                    Vec3 f_ab = eval_deon(wi, wo, t, p);
                    Vec3 f_ba = eval_deon(wo, wi, t, p);
                    float err = std::max(
                        {fabsf(f_ab.x - f_ba.x), fabsf(f_ab.y - f_ba.y), fabsf(f_ab.z - f_ba.z)});
                    if (err > 1e-3F) {
                        spdlog::error(
                            "  FAIL reciprocity err={:.5f}  ti={:.2f} pi={:.2f} to={:.2f} po={:.2f}",
                            err,
                            ti,
                            pi,
                            to,
                            po);
                        passed = false;
                    }
                }
}
}
}
    return passed;
}

static bool test_r_lobe_specular() {
    DeonParams p;
    p.sigma_a    = {0.F, 0.F, 0.F};
    p.beta_TT    = 1e-6F;
    p.beta_TRT   = 1e-6F;
    const Vec3 t = {0.F, 0.F, 1.F};

    Vec3 wi      = dir_from_angles(0.3F, 0.F);
    Vec3 wo_spec = dir_from_angles(-0.3F, 0.F);
    Vec3 wo_off  = dir_from_angles(0.3F, 2.F);

    Vec3 f_spec = eval_deon(wi, wo_spec, t, p);
    Vec3 f_off  = eval_deon(wi, wo_off, t, p);

    bool ok = (f_spec.x > f_off.x) && (f_spec.y > f_off.y) && (f_spec.z > f_off.z);
    if (!ok) {
        spdlog::error("  FAIL R_lobe: f_spec=({:.4f},{:.4f},{:.4f}) f_off=({:.4f},{:.4f},{:.4f})",
                      f_spec.x,
                      f_spec.y,
                      f_spec.z,
                      f_off.x,
                      f_off.y,
                      f_off.z);
}
    return ok;
}

int main() {
    spdlog::set_pattern("[%T.%e] [%^%l%$] %v");
    spdlog::info("=== d'Eon BSDF unit tests ===");

    bool ok1 = test_white_furnace();
    spdlog::log(ok1 ? spdlog::level::info : spdlog::level::err,
                "[{}] white_furnace",
                ok1 ? "PASS" : "FAIL");

    bool ok2 = test_reciprocity();
    spdlog::log(
        ok2 ? spdlog::level::info : spdlog::level::err, "[{}] reciprocity", ok2 ? "PASS" : "FAIL");

    bool ok3 = test_r_lobe_specular();
    spdlog::log(ok3 ? spdlog::level::info : spdlog::level::err,
                "[{}] R_lobe_specular",
                ok3 ? "PASS" : "FAIL");

    bool all = ok1 && ok2 && ok3;
    spdlog::log(all ? spdlog::level::info : spdlog::level::err,
                all ? "All tests PASSED." : "Some tests FAILED.");
    return all ? 0 : 1;
}
