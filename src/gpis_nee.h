#pragma once
#include "math.h"

namespace m3hair {

// ---------------------------------------------------------------------------
// Hit information needed for GPIS next-event estimation (§6, Xu et al. 2025)
// ---------------------------------------------------------------------------
struct GpisHitInfo {
    vec3 p;       // hit point x_s
    vec3 grad_mu; // ∇μ(x_s)  — mean-field gradient
    float gz;     // directional noise derivative ∂ψ/∂z at x_s (ray dir z)
    float l2;     // l² = (cell_size/3)²
    float kdd0;   // κ''(0) = −1/(2l²)
};

// ---------------------------------------------------------------------------
// Gradient-distribution helpers (Eq. 24, isotropic stationary case)
// ---------------------------------------------------------------------------

// PDF of the 2D transverse gradient under the isotropic Gaussian prior.
//   p(gxy) = (kdd0_abs / π) · exp(−kdd0_abs · ‖gxy‖²)
//   kdd0_abs = |κ''(0)| = 1/(2l²)
float eval_gxy_pdf(vec2 gxy, float kdd0_abs);

// Sample gxy ~ N(0, 1/(2l²) · I₂) via Box-Muller.
//   kdd0_abs = 1/(2l²);  u1, u2 ∈ (0,1)
vec2 sample_gxy(float kdd0_abs, float u1, float u2);

// Convert a desired outgoing direction omega_s and the known gz component to
// the implied transverse gradient gxy (Fig. 11 of the paper).
// n = normalize(grad_mu + (gxy, gz)); given omega_s we want n = n_target.
// Returns the gxy that places the normal along the half-vector direction.
vec2 to_gxy(vec3 grad_mu, vec3 omega_s, float gz);

// Convert a solid-angle light PDF p_omega to the gradient-space measure p_g.
//   p^g(gxy) = p^n(n) · |n.z| / ‖g‖²,    g = (gxy, gz)
//   p^n(n)   = 4·(n·omega_s)·p_omega
// Returns 0 if the measure conversion is degenerate.
float to_gxy_measure(float p_omega, vec3 omega_s, vec3 g);

// Balance-heuristic MIS weight for two strategies with PDFs p_a and p_b.
float mis_weight(float p_a, float p_b);

// ---------------------------------------------------------------------------
// Full NEE contribution at a GPIS shading point (Algorithm 4, §6)
//
//  hit      — shading-point data (from extract_hit_info)
//  omega_o  — outgoing direction toward camera (unit)
//  omega_s  — sampled light direction (unit)
//  L_e      — emitter radiance (solid-angle)
//  p_light  — solid-angle PDF of omega_s
//  rho      — BRDF × cosine factor evaluated by the caller
//  u_gxy1,2 — uniform samples in (0,1) for the gradient prior sampling
//
// Returns the MIS-weighted radiance contribution (multiply by visibility).
// ---------------------------------------------------------------------------
vec3 eval_nee_contribution(const GpisHitInfo &hit,
                           vec3 omega_o,
                           vec3 omega_s,
                           vec3 L_e,
                           float p_light,
                           float rho,
                           float u_gxy1,
                           float u_gxy2);

} // namespace m3hair
