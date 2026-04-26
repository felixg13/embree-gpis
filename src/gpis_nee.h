#pragma once
#include "math.h"

namespace m3hair {

struct GpisHitInfo {
    Vec3 p;
    Vec3 grad_mu;
    float gz;
    float l2;
    float kdd0;
};

float eval_gxy_pdf(Vec2 gxy, float kdd0_abs);
Vec2 sample_gxy(float kdd0_abs, float u1, float u2);
Vec2 to_gxy(Vec3 grad_mu, Vec3 omega_s, float gz);
float to_gxy_measure(float p_omega, Vec3 omega_s, Vec3 g);
float mis_weight(float p_a, float p_b);
Vec3 eval_nee_contribution(const GpisHitInfo &hit,
                           Vec3 omega_o,
                           Vec3 omega_s,
                           Vec3 L_e,
                           float p_light,
                           float rho,
                           float u_gxy1 = 0.F,
                           float u_gxy2 = 0.F);

} // namespace m3hair
