// GPIS Next-Event Estimation helpers — Xu et al. 2025, §6 / Algorithm 4.
// Implements gradient-distribution importance sampling for GPIS surfaces.

#include "gpis_nee.h"

#include <algorithm>
#include <cmath>

namespace m3hair {

// ---------------------------------------------------------------------------
// eval_gxy_pdf
// p(gxy) = (kdd0_abs / π) · exp(−kdd0_abs · ‖gxy‖²)
// This is a 2D isotropic Gaussian with variance σ² = 1/(2·kdd0_abs).
// ---------------------------------------------------------------------------
float eval_gxy_pdf(vec2 gxy, float kdd0_abs) {
    float g2 = dot(gxy, gxy);
    return (kdd0_abs / 3.14159265f) * expf(-kdd0_abs * g2);
}

// ---------------------------------------------------------------------------
// sample_gxy — Box-Muller from N(0, 1/(2·kdd0_abs) · I₂)
// σ² = 1/(2·kdd0_abs)  →  σ = 1/sqrt(2·kdd0_abs)
// ---------------------------------------------------------------------------
vec2 sample_gxy(float kdd0_abs, float u1, float u2) {
    // Clamp u1 away from 0 to avoid log(0)
    u1                = std::max(u1, 1e-7f);
    const float sigma = 1.f / sqrtf(2.f * kdd0_abs);
    const float r     = sigma * sqrtf(-2.f * logf(u1));
    const float theta = 2.f * 3.14159265f * u2;
    return {r * cosf(theta), r * sinf(theta)};
}

// ---------------------------------------------------------------------------
// to_gxy — find the transverse gradient (gxy) that places the surface normal
// along the half-vector n_target = normalize(omega_o + omega_s).
//
// The surface normal is n = normalize(grad_mu + g),  g = (gxy.x, gxy.y, gz).
// We want n = n_target.  Since n must be proportional to grad_mu + g:
//   grad_mu + (gxy, gz) = s · n_target   for some s > 0
//   (gxy, gz) = s · n_target - grad_mu
// The z-component constraint fixes s:
//   gz = s · n_target.z - grad_mu.z   →   s = (gz + grad_mu.z) / n_target.z
//
// If n_target.z ≈ 0 the inversion is degenerate; return (0,0) in that case.
// ---------------------------------------------------------------------------
vec2 to_gxy(vec3 grad_mu, vec3 omega_s, float gz) {
    // Half-vector (un-normalized) used as the target normal direction
    // In a mirror-BRDF context, n_target = normalize(omega_o + omega_s).
    // Here omega_s is the light direction and omega_o is baked into the call —
    // the caller passes omega_s directly; we treat it as the target normal dir.
    vec3 n_target = normalize(omega_s);

    if (fabsf(n_target.z) < 1e-6f)
        return {0.f, 0.f};

    float s  = (gz + grad_mu.z) / n_target.z;
    float gx = s * n_target.x - grad_mu.x;
    float gy = s * n_target.y - grad_mu.y;
    return {gx, gy};
}

// ---------------------------------------------------------------------------
// to_gxy_measure — change of variables from solid-angle PDF to gradient PDF.
//
// The normal is n = normalize(grad_mu + g),  g = (gxy, gz).
// The Jacobian from solid-angle to gradient measure is derived in §6:
//   p^n(n) = 4·(n·omega_s)·p^omega(omega_s)
//   p^g(g) = p^n(n) · |∂n/∂g|
// For an isotropic surface, |∂n/∂g| = 1/‖g‖²  (up to the |n_z| projection).
// Returns 0 for degenerate configurations (g ≈ 0 or backfacing).
// ---------------------------------------------------------------------------
float to_gxy_measure(float p_omega, vec3 omega_s, vec3 g) {
    float g_len2 = dot(g, g);
    if (g_len2 < 1e-12f)
        return 0.f;

    vec3 n        = normalize(g); // g = grad_mu + delta, effectively n direction
    float n_dot_s = dot(n, omega_s);
    if (n_dot_s <= 0.f)
        return 0.f;

    float p_n = 4.f * n_dot_s * p_omega;
    return p_n * fabsf(n.z) / g_len2;
}

// ---------------------------------------------------------------------------
// mis_weight — balance heuristic for two strategies
// ---------------------------------------------------------------------------
float mis_weight(float p_a, float p_b) {
    float sum = p_a + p_b;
    if (sum < 1e-30f)
        return 0.f;
    return p_a / sum;
}

// ---------------------------------------------------------------------------
// eval_nee_contribution — Algorithm 4 (§6, Xu et al. 2025)
// ---------------------------------------------------------------------------
vec3 eval_nee_contribution(const GpisHitInfo &hit,
                           vec3 omega_o,
                           vec3 omega_s,
                           vec3 L_e,
                           float p_light,
                           float rho,
                           float u_gxy1,
                           float u_gxy2) {
    (void)omega_o; // used by caller for BRDF; rho already incorporates it
    (void)u_gxy1;
    (void)u_gxy2;

    const float kdd0_abs = -hit.kdd0; // |κ''(0)| = 1/(2l²), positive

    // 1. Known gz component (ray-direction projection of noise gradient)
    float gz = hit.gz;

    // 2. Implied transverse gradient for the NEE direction
    vec2 gxy_nee = to_gxy(hit.grad_mu, omega_s, gz);
    vec3 g_nee   = {gxy_nee.x, gxy_nee.y, gz};

    // 3. PDF of this g_nee under the gradient-space light measure
    float p_nee_gxy = to_gxy_measure(p_light, omega_s, g_nee);

    // 4. PDF of g_nee under the gradient prior
    float p_gp_gxy = eval_gxy_pdf(gxy_nee, kdd0_abs);

    // 5. Surface normal implied by g_nee
    //    (used for shading; rho already computed by caller)
    //    vec3 n_nee = normalize(hit.grad_mu + g_nee);

    // 6. MIS weight (balance heuristic, NEE vs. gradient-prior sampling)
    float w = mis_weight(p_nee_gxy, p_gp_gxy);

    // 7. Guard against degenerate measure
    float denom = p_nee_gxy + p_gp_gxy;
    if (denom < 1e-30f)
        return {0.f, 0.f, 0.f};

    // Contribution = L_e · rho · w / denom  (visibility applied by caller)
    float weight = rho * w / denom;
    return L_e * weight;
}

} // namespace m3hair
