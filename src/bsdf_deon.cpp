#include "bsdf_deon.h"

#include <algorithm>
#include <cmath>

namespace m3hair {

static constexpr float PI      = 3.14159265358979f;
static constexpr float SQRT2PI = 2.50662827463f;

static float fresnel(float eta, float cos_i) {
    cos_i        = std::abs(cos_i);
    float sin_i2 = std::max(0.f, 1.f - cos_i * cos_i);
    float sin_t2 = sin_i2 / (eta * eta);
    if (sin_t2 >= 1.f)
        return 1.f;
    float cos_t = sqrtf(1.f - sin_t2);
    float rs    = (cos_i - eta * cos_t) / (cos_i + eta * cos_t);
    float rp    = (eta * cos_i - cos_t) / (eta * cos_i + cos_t);
    return 0.5f * (rs * rs + rp * rp);
}

static float gauss(float beta, float x) {
    return expf(-x * x / (2.f * beta * beta)) / (beta * SQRT2PI);
}

static float wrap(float a) {
    a = fmodf(a, 2.f * PI);
    if (a > PI)
        a -= 2.f * PI;
    if (a < -PI)
        a += 2.f * PI;
    return a;
}

// azimuthal exit angle for lobe p, etap = sqrt(eta^2 - sin^2(theta_d)) / cos(theta_d)
static float phi_exit(int p, float h, float etap) {
    float gamma_i = asinf(std::clamp(h, -1.f, 1.f));
    float gamma_t = asinf(std::clamp(h / etap, -1.f, 1.f));
    return 2.f * p * gamma_t - 2.f * gamma_i + (float)p * PI;
}

// 16-point midpoint quadrature over fiber cross-section
static vec3 eval_Np(int p, float phi, float etap, float cos_theta_d, const vec3 &sa, float beta_n) {
    const int N = 16;
    vec3 acc(0.f);
    for (int j = 0; j < N; ++j) {
        float h      = -1.f + (2.f * j + 1.f) / N;
        float cos_gi = sqrtf(std::max(0.f, 1.f - h * h));
        float F      = fresnel(etap, cos_gi);
        // absorption over chord through fiber
        vec3 T = {expf(-2.f * sa.x * cos_gi / cos_theta_d),
                  expf(-2.f * sa.y * cos_gi / cos_theta_d),
                  expf(-2.f * sa.z * cos_gi / cos_theta_d)};
        vec3 A;
        if (p == 0) {
            A = vec3(F);
        } else if (p == 1) {
            float fac = (1.f - F) * (1.f - F);
            A         = T * fac;
        } else {
            float fac = (1.f - F) * (1.f - F) * F;
            A         = T * T * fac;
        }
        float g = gauss(beta_n, wrap(phi - phi_exit(p, h, etap)));
        acc     = acc + A * g;
    }
    return acc * (1.f / N);
}

vec3 eval_deon(const vec3 &wi, const vec3 &wo, const vec3 &T, const DeonParams &p) {
    float si = std::clamp(dot(wi, T), -1.f, 1.f);
    float so = std::clamp(dot(wo, T), -1.f, 1.f);
    float ti = asinf(si), to = asinf(so);
    float td = (to - ti) * 0.5f, th = (to + ti) * 0.5f;
    float cos_td  = cosf(td);
    float cos2_td = cos_td * cos_td;
    if (cos2_td < 1e-7f)
        return vec3(0.f);

    vec3 wi_p = wi - si * T, wo_p = wo - so * T;
    float li = length(wi_p), lo = length(wo_p);
    float phi = 0.f;
    if (li > 1e-5f && lo > 1e-5f)
        phi = acosf(std::clamp(dot(wi_p, wo_p) / (li * lo), -1.f, 1.f));

    float sin_td2 = 1.f - cos2_td;
    float etap    = sqrtf(std::max(0.f, p.eta * p.eta - sin_td2)) / cos_td;
    etap          = std::max(etap, 1.001f);

    const float alpha[3]  = {-p.alpha, 0.5f * p.alpha, -1.5f * p.alpha};
    const float beta_p[3] = {p.beta_R, p.beta_TT, p.beta_TRT};

    vec3 result(0.f);
    for (int lobe = 0; lobe < 3; ++lobe) {
        // 0.5 corrects for M being normalized in theta_h while integral runs over theta_o
        float M = 0.5f * gauss(beta_p[lobe], th - alpha[lobe]);
        vec3 N  = eval_Np(lobe, phi, etap, cos_td, p.sigma_a, p.beta_n);
        result  = result + N * M;
    }
    return result / cos2_td;
}

static float lum(const vec3 &v) {
    return 0.2126f * v.x + 0.7152f * v.y + 0.0722f * v.z;
}

vec3 sample_deon(const vec3 &wo,
                 const vec3 &T,
                 const DeonParams &p,
                 float u1,
                 float u2,
                 float u3,
                 float u4,
                 vec3 &wi_out,
                 float &pdf_out) {
    float so = std::clamp(dot(wo, T), -1.f, 1.f);
    float to = asinf(so);

    float F0   = fresnel(p.eta, 1.f);
    float w[3] = {F0, (1.f - F0) * (1.f - F0), (1.f - F0) * (1.f - F0) * F0};
    float wsum = w[0] + w[1] + w[2];
    if (wsum < 1e-7f) {
        pdf_out = 0.f;
        return vec3(0.f);
    }
    for (float &x : w)
        x /= wsum;

    int lobe  = 2;
    float cdf = 0.f;
    for (int i = 0; i < 3; ++i) {
        cdf += w[i];
        if (u1 < cdf) {
            lobe = i;
            break;
        }
    }

    const float alpha_p[3] = {-p.alpha, 0.5f * p.alpha, -1.5f * p.alpha};
    const float beta_p[3]  = {p.beta_R, p.beta_TT, p.beta_TRT};

    u2       = std::max(u2, 1e-7f);
    float z  = sqrtf(-2.f * logf(u2)) * cosf(2.f * PI * u3);
    float th = alpha_p[lobe] + beta_p[lobe] * z;
    float ti = std::clamp(2.f * th - to, -PI * 0.5f, PI * 0.5f);
    float si = sinf(ti);
    float ci = sqrtf(std::max(0.f, 1.f - si * si));

    float td     = (to - ti) * 0.5f;
    float cos_td = cosf(td);
    if (cos_td < 1e-5f) {
        pdf_out = 0.f;
        return vec3(0.f);
    }

    float etap = sqrtf(std::max(0.f, p.eta * p.eta - (1.f - cos_td * cos_td))) / cos_td;
    etap       = std::max(etap, 1.001f);

    const int NPHI = 64;
    float table[NPHI];
    float tsum = 0.f;
    for (int k = 0; k < NPHI; ++k) {
        float phi_k = (k + 0.5f) * (2.f * PI / NPHI);
        table[k]    = std::max(0.f, lum(eval_Np(lobe, phi_k, etap, cos_td, p.sigma_a, p.beta_n)));
        tsum += table[k];
    }
    float phi_s = 0.f;
    if (tsum < 1e-7f) {
        phi_s = u4 * 2.f * PI;
    } else {
        float target = u4 * tsum, acc = 0.f;
        phi_s = (NPHI - 0.5f) * (2.f * PI / NPHI);
        for (int k = 0; k < NPHI; ++k) {
            acc += table[k];
            if (acc >= target) {
                phi_s = (k + 0.5f) * (2.f * PI / NPHI);
                break;
            }
        }
    }

    vec3 wo_p = wo - so * T;
    float lop = length(wo_p);
    vec3 phi_x;
    if (lop > 1e-5f) {
        phi_x = wo_p / lop;
    } else {
        vec3 arb = (fabsf(T.x) < 0.9f) ? vec3(1, 0, 0) : vec3(0, 1, 0);
        phi_x    = normalize(cross(T, arb));
    }
    vec3 phi_y = cross(T, phi_x);

    wi_out = normalize(T * si + phi_x * (ci * cosf(phi_s)) + phi_y * (ci * sinf(phi_s)));

    float M_pdf  = gauss(beta_p[lobe], th - alpha_p[lobe]);
    float Np_val = (tsum > 1e-7f) ? lum(eval_Np(lobe, phi_s, etap, cos_td, p.sigma_a, p.beta_n)) /
                                        (tsum * (2.f * PI / NPHI))
                                  : 1.f / (2.f * PI);

    pdf_out = w[lobe] * M_pdf * 0.5f * Np_val;
    if (pdf_out < 1e-7f) {
        pdf_out = 0.f;
        return vec3(0.f);
    }

    return eval_deon(wi_out, wo, T, p);
}

} // namespace m3hair
