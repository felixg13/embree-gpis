#include "bsdf_deon.h"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace m3hair {

static constexpr float PI      = std::numbers::pi_v<float>;
static constexpr float SQRT2PI = 2.50662827463F;

static float fresnel(float eta, float cos_i) {
    cos_i        = std::abs(cos_i);
    float sin_i2 = std::max(0.F, 1.F - (cos_i * cos_i));
    float sin_t2 = sin_i2 / (eta * eta);
    if (sin_t2 >= 1.F) {
        return 1.F;
}
    float cos_t = sqrtf(1.F - sin_t2);
    float rs    = (cos_i - (eta * cos_t)) / (cos_i + (eta * cos_t));
    float rp    = ((eta * cos_i) - cos_t) / ((eta * cos_i) + cos_t);
    return 0.5F * ((rs * rs) + (rp * rp));
}

static float gauss(float beta, float x) {
    return expf(-x * x / (2.F * beta * beta)) / (beta * SQRT2PI);
}

static float wrap(float a) {
    a = fmodf(a, 2.F * PI);
    if (a > PI) {
        a -= 2.F * PI;
}
    if (a < -PI) {
        a += 2.F * PI;
}
    return a;
}

// azimuthal exit angle for lobe p, etap = sqrt(eta^2 - sin^2(theta_d)) / cos(theta_d)
static float phi_exit(int p, float h, float etap) {
    float gamma_i = asinf(std::clamp(h, -1.F, 1.F));
    float gamma_t = asinf(std::clamp(h / etap, -1.F, 1.F));
    return (2.F * p * gamma_t) - (2.F * gamma_i) + (static_cast<float>(p) * PI);
}

// 16-point midpoint quadrature over fiber cross-section
static Vec3 eval_np(int p, float phi, float etap, float cos_theta_d, const Vec3 &sa, float beta_n) {
    const int n = 16;
    Vec3 acc(0.F);
    for (int j = 0; j < n; ++j) {
        float h      = -1.F + (((2.F * j) + 1.F) / n);
        float cos_gi = sqrtf(std::max(0.F, 1.F - (h * h)));
        float f      = fresnel(etap, cos_gi);
        // absorption over chord through fiber
        Vec3 t = {expf(-2.F * sa.x * cos_gi / cos_theta_d),
                  expf(-2.F * sa.y * cos_gi / cos_theta_d),
                  expf(-2.F * sa.z * cos_gi / cos_theta_d)};
        Vec3 a;
        if (p == 0) {
            a = Vec3(f);
        } else if (p == 1) {
            float fac = (1.F - f) * (1.F - f);
            a         = t * fac;
        } else {
            float fac = (1.F - f) * (1.F - f) * f;
            a         = t * t * fac;
        }
        float g = gauss(beta_n, wrap(phi - phi_exit(p, h, etap)));
        acc     = acc + a * g;
    }
    return acc * (1.F / n);
}

Vec3 eval_deon(const Vec3 &wi, const Vec3 &wo, const Vec3 &T, const DeonParams &p) {
    float si = std::clamp(dot(wi, T), -1.F, 1.F);
    float so = std::clamp(dot(wo, T), -1.F, 1.F);
    float ti = asinf(si);
    float to = asinf(so);
    float td = (to - ti) * 0.5f;
    float th = (to + ti) * 0.5f;
    float cos_td  = cosf(td);
    float cos2_td = cos_td * cos_td;
    if (cos2_td < 1e-7F) {
        return Vec3(0.F);
}

    Vec3 wi_p = wi - si * T;
    Vec3 wo_p = wo - so * T;
    float li = length(wi_p);
    float lo = length(wo_p);
    float phi = 0.F;
    if (li > 1e-5F && lo > 1e-5F) {
        phi = acosf(std::clamp(dot(wi_p, wo_p) / (li * lo), -1.F, 1.F));
}

    float sin_td2 = 1.F - cos2_td;
    float etap    = sqrtf(std::max(0.F, (p.eta * p.eta) - sin_td2)) / cos_td;
    etap          = std::max(etap, 1.001F);

    const float alpha[3]  = {-p.alpha, 0.5F * p.alpha, -1.5F * p.alpha};
    const float beta_p[3] = {p.beta_R, p.beta_TT, p.beta_TRT};

    Vec3 result(0.F);
    for (int lobe = 0; lobe < 3; ++lobe) {
        // 0.5 corrects for M being normalized in theta_h while integral runs over theta_o
        float m = 0.5F * gauss(beta_p[lobe], th - alpha[lobe]);
        Vec3 n  = eval_np(lobe, phi, etap, cos_td, p.sigma_a, p.beta_n);
        result  = result + n * m;
    }
    return result / cos2_td;
}

static float lum(const Vec3 &v) {
    return (0.2126F * v.x) + (0.7152F * v.y) + (0.0722F * v.z);
}

Vec3 sample_deon(const Vec3 &wo,
                 const Vec3 &T,
                 const DeonParams &p,
                 float u1,
                 float u2,
                 float u3,
                 float u4,
                 Vec3 &wi_out,
                 float &pdf_out) {
    float so = std::clamp(dot(wo, T), -1.F, 1.F);
    float to = asinf(so);

    float f0   = fresnel(p.eta, 1.F);
    float w[3] = {f0, (1.F - f0) * (1.F - f0), (1.F - f0) * (1.F - f0) * f0};
    float wsum = w[0] + w[1] + w[2];
    if (wsum < 1e-7F) {
        pdf_out = 0.F;
        return Vec3(0.F);
    }
    for (float &x : w) {
        x /= wsum;
}

    int lobe  = 2;
    float cdf = 0.F;
    for (int i = 0; i < 3; ++i) {
        cdf += w[i];
        if (u1 < cdf) {
            lobe = i;
            break;
        }
    }

    const float alpha_p[3] = {-p.alpha, 0.5F * p.alpha, -1.5F * p.alpha};
    const float beta_p[3]  = {p.beta_R, p.beta_TT, p.beta_TRT};

    u2       = std::max(u2, 1e-7F);
    float z  = sqrtf(-2.F * logf(u2)) * cosf(2.F * PI * u3);
    float th = alpha_p[lobe] + (beta_p[lobe] * z);
    float ti = std::clamp((2.F * th) - to, -PI * 0.5F, PI * 0.5F);
    float si = sinf(ti);
    float ci = sqrtf(std::max(0.F, 1.F - (si * si)));

    float td     = (to - ti) * 0.5F;
    float cos_td = cosf(td);
    if (cos_td < 1e-5F) {
        pdf_out = 0.F;
        return Vec3(0.F);
    }

    float etap = sqrtf(std::max(0.F, (p.eta * p.eta) - (1.F - (cos_td * cos_td)))) / cos_td;
    etap       = std::max(etap, 1.001F);

    const int nphi = 64;
    float table[nphi];
    float tsum = 0.F;
    for (int k = 0; k < nphi; ++k) {
        float phi_k = (k + 0.5F) * (2.F * PI / nphi);
        table[k]    = std::max(0.F, lum(eval_np(lobe, phi_k, etap, cos_td, p.sigma_a, p.beta_n)));
        tsum += table[k];
    }
    float phi_s = 0.F;
    if (tsum < 1e-7F) {
        phi_s = u4 * 2.F * PI;
    } else {
        float target = u4 * tsum;
        float acc = 0.f;
        phi_s = (nphi - 0.5F) * (2.F * PI / nphi);
        for (int k = 0; k < nphi; ++k) {
            acc += table[k];
            if (acc >= target) {
                phi_s = (k + 0.5F) * (2.F * PI / nphi);
                break;
            }
        }
    }

    Vec3 wo_p = wo - so * T;
    float lop = length(wo_p);
    Vec3 phi_x;
    if (lop > 1e-5F) {
        phi_x = wo_p / lop;
    } else {
        Vec3 arb = (fabsf(T.x) < 0.9F) ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
        phi_x    = normalize(cross(T, arb));
    }
    Vec3 phi_y = cross(T, phi_x);

    wi_out = normalize(T * si + phi_x * (ci * cosf(phi_s)) + phi_y * (ci * sinf(phi_s)));

    float m_pdf  = gauss(beta_p[lobe], th - alpha_p[lobe]);
    float np_val = (tsum > 1e-7F) ? lum(eval_np(lobe, phi_s, etap, cos_td, p.sigma_a, p.beta_n)) /
                                        (tsum * (2.F * PI / nphi))
                                  : 1.F / (2.F * PI);

    pdf_out = w[lobe] * m_pdf * 0.5F * np_val;
    if (pdf_out < 1e-7F) {
        pdf_out = 0.F;
        return Vec3(0.F);
    }

    return eval_deon(wi_out, wo, T, p);
}

} // namespace m3hair
