#pragma once
#include "math.h"

namespace m3hair {

struct DeonParams {
    float beta_R   = 0.10f;
    float beta_TT  = 0.05f;
    float beta_TRT = 0.20f;
    float alpha    = 0.035f;
    float eta      = 1.55f;
    vec3 sigma_a   = {0.25f, 0.7f, 2.0f};
    float beta_n   = 0.30f;
};

vec3 eval_deon(const vec3 &wi, const vec3 &wo, const vec3 &T, const DeonParams &p);

vec3 sample_deon(const vec3 &wo,
                 const vec3 &T,
                 const DeonParams &p,
                 float u1,
                 float u2,
                 float u3,
                 float u4,
                 vec3 &wi_out,
                 float &pdf_out);

} // namespace m3hair
