#pragma once
#include "math.h"

namespace m3hair {

struct GpisHitInfo {
    vec3 p;
    vec3 grad_mu;
    float gz;
    float l2;
    float kdd0;
};

float eval_gxy_pdf(vec2 gxy, float kdd0_abs);
vec2 sample_gxy(float kdd0_abs, float u1, float u2);
vec2 to_gxy(vec3 grad_mu, vec3 omega_s, float gz);
float to_gxy_measure(float p_omega, vec3 omega_s, vec3 g);
float mis_weight(float p_a, float p_b);
vec3 eval_nee_contribution(const GpisHitInfo &hit,
                           vec3 omega_o,
                           vec3 omega_s,
                           vec3 L_e,
                           float p_light,
                           float rho,
                           float u_gxy1 = 0.f,
                           float u_gxy2 = 0.f);

} // namespace m3hair
