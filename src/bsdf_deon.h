#pragma once
#include "math.h"

namespace m3hair {

struct DeonParams {
    float beta_R   = 0.10f;   // longitudinal roughness R
    float beta_TT  = 0.05f;   // = beta_R / 2
    float beta_TRT = 0.20f;   // = 2 * beta_R
    float alpha    = 0.035f;  // tilt angle (radians, ~2 deg)
    float eta      = 1.55f;   // IOR
    vec3  sigma_a  = {0.25f, 0.7f, 2.0f};  // golden default
    float beta_n   = 0.30f;   // azimuthal roughness
};

// Evaluate f(wi, wo) given fiber tangent T.
// Both wi and wo point away from the surface (outward convention).
vec3 eval_deon(const vec3& wi, const vec3& wo, const vec3& T,
               const DeonParams& p);

// Importance-sample wi given wo.
// u1        : lobe selection
// u2, u3    : Box-Muller pair for theta_i
// u4        : phi CDF lookup
// Returns f(wi_out, wo); sets wi_out and pdf_out.
vec3 sample_deon(const vec3& wo, const vec3& T, const DeonParams& p,
                 float u1, float u2, float u3, float u4,
                 vec3& wi_out, float& pdf_out);

} // namespace m3hair
