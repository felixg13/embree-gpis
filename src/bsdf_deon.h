#pragma once
#include "math.h"

namespace m3hair {

struct DeonParams {
    float beta_R   = 0.10F;
    float beta_TT  = 0.05F;
    float beta_TRT = 0.20F;
    float alpha    = 0.035F;
    float eta      = 1.55F;
    Vec3 sigma_a   = {0.25F, 0.7F, 2.0F};
    float beta_n   = 0.30F;
};

Vec3 eval_deon(const Vec3 &wi, const Vec3 &wo, const Vec3 &T, const DeonParams &p);

Vec3 sample_deon(const Vec3 &wo,
                 const Vec3 &T,
                 const DeonParams &p,
                 float u1,
                 float u2,
                 float u3,
                 float u4,
                 Vec3 &wi_out,
                 float &pdf_out);

} // namespace m3hair
