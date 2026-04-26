#include "gpis_nee.h"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace m3hair {

float eval_gxy_pdf(Vec2 gxy, float kdd0_abs) {
    float g2 = dot(gxy, gxy);
    return (kdd0_abs / std::numbers::pi_v<float>) * expf(-kdd0_abs * g2);
}

Vec2 sample_gxy(float kdd0_abs, float u1, float u2) {
    u1                = std::max(u1, 1e-7F);
    const float sigma = 1.F / sqrtf(2.F * kdd0_abs);
    const float r     = sigma * sqrtf(-2.F * logf(u1));
    const float theta = 2.F * std::numbers::pi_v<float> * u2;
    return {r * cosf(theta), r * sinf(theta)};
}

// find gxy such that normalize(grad_mu + (gxy, gz)) = n_target
Vec2 to_gxy(Vec3 grad_mu, Vec3 omega_s, float gz) {
    Vec3 n_target = normalize(omega_s);

    if (fabsf(n_target.z) < 1e-6F) {
        return {0.F, 0.F};
}

    float s  = (gz + grad_mu.z) / n_target.z;
    float gx = (s * n_target.x) - grad_mu.x;
    float gy = (s * n_target.y) - grad_mu.y;
    return {gx, gy};
}

float to_gxy_measure(float p_omega, Vec3 omega_s, Vec3 g) {
    float g_len2 = dot(g, g);
    if (g_len2 < 1e-12F) {
        return 0.F;
}

    Vec3 n        = normalize(g);
    float n_dot_s = dot(n, omega_s);
    if (n_dot_s <= 0.F) {
        return 0.F;
}

    float p_n = 4.F * n_dot_s * p_omega;
    return p_n * fabsf(n.z) / g_len2;
}

float mis_weight(float p_a, float p_b) {
    float sum = p_a + p_b;
    if (sum < 1e-30F) {
        return 0.F;
}
    return p_a / sum;
}

Vec3 eval_nee_contribution(const GpisHitInfo &hit,
                           Vec3 /*omega_o*/,
                           Vec3 omega_s,
                           Vec3 L_e,
                           float p_light,
                           float rho,
                           float /*u_gxy1*/,
                           float /*u_gxy2*/) {

    const float kdd0_abs = -hit.kdd0;

    float gz         = hit.gz;
    Vec2 gxy_nee     = to_gxy(hit.grad_mu, omega_s, gz);
    Vec3 g_nee       = {gxy_nee.x, gxy_nee.y, gz};
    float p_nee_gxy  = to_gxy_measure(p_light, omega_s, g_nee);
    float p_gp_gxy   = eval_gxy_pdf(gxy_nee, kdd0_abs);
    float w          = mis_weight(p_nee_gxy, p_gp_gxy);
    float denom      = p_nee_gxy + p_gp_gxy;
    if (denom < 1e-30F) {
        return {0.F, 0.F, 0.F};
}

    float weight = rho * w / denom;
    return L_e * weight;
}

} // namespace m3hair
