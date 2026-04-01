// GPIS (Gaussian Process Implicit Surface) user geometry for Embree 4.
// Implements sparse convolution noise following Xu et al. 2025.
//
// The noise field is:
//   f(x) = μ(x) + A · scnoise(x)
//
// where scnoise is a sparse convolution noise (Lagae et al. 2011):
//   scnoise(x) = (1/√Var) · Σ_cells Σ_k  w_k · φ(x − x_k)
//
//   • x_k: impulse positions, uniformly distributed in each grid cell
//   • w_k: ±1 Bernoulli weights
//   • φ: Gaussian kernel exp(−‖x−x_k‖²/(2l²)), l = cell_size/3
//   • Var = N · π^(3/2) · l³ / l³  (analytic for Gaussian kernel, see below)
//
// The gradient ∇f is computed analytically:
//   ∇f(x) = ∇μ(x) + A · ∇scnoise(x)
//   ∇φ(x−x_k) = −(x−x_k)/l² · φ(x−x_k)
//
// Each cell's impulse positions and weights are seeded deterministically from
// the cell's integer coordinates + a global seed, giving a continuous coherent
// noise field regardless of which primitive or ray evaluates it.

#include "user_geo.h"

#include "gpis_nee.h"
#include "math.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace m3hair {

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

static vec3 from_f4(const float4 &v) {
    return {v.x, v.y, v.z};
}

static vec3 bspline_pos(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float t) {
    float t2 = t * t, t3 = t2 * t;
    float b0 = (-t3 + 3 * t2 - 3 * t + 1) * (1.f / 6.f);
    float b1 = (3 * t3 - 6 * t2 + 4) * (1.f / 6.f);
    float b2 = (-3 * t3 + 3 * t2 + 3 * t + 1) * (1.f / 6.f);
    float b3 = t3 * (1.f / 6.f);
    return p0 * b0 + p1 * b1 + p2 * b2 + p3 * b3;
}

// Approximate closest t on B-spline to point x via 12-sample search
static float closest_t(vec3 p0, vec3 p1, vec3 p2, vec3 p3, vec3 x) {
    float best_t = 0.f, best_d2 = 1e30f;
    for (int i = 0; i <= 12; ++i) {
        float t  = i / 12.f;
        vec3 pt  = bspline_pos(p0, p1, p2, p3, t);
        float d2 = dot(pt - x, pt - x);
        if (d2 < best_d2) {
            best_d2 = d2;
            best_t  = t;
        }
    }
    float eps = 1.f / 24.f;
    for (int iter = 0; iter < 2; ++iter) {
        vec3 pa  = bspline_pos(p0, p1, p2, p3, std::max(0.f, best_t - eps));
        vec3 pb  = bspline_pos(p0, p1, p2, p3, std::min(1.f, best_t + eps));
        float da = dot(pa - x, pa - x), db = dot(pb - x, pb - x);
        if (da < best_d2) {
            best_d2 = da;
            best_t  = std::max(0.f, best_t - eps);
        }
        if (db < best_d2) {
            best_d2 = db;
            best_t  = std::min(1.f, best_t + eps);
        }
        eps *= 0.5f;
    }
    return best_t;
}

// Mean field value only (cheap path, used during marching)
static float mean_field_val(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float r, vec3 x) {
    float t = closest_t(p0, p1, p2, p3, x);
    vec3 cp = bspline_pos(p0, p1, p2, p3, t);
    vec3 d  = x - cp;
    return r - sqrtf(dot(d, d));
}

// Mean field value + gradient (used only at the hit point for normal)
static float mean_field(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float r, vec3 x, vec3 &out_grad_mu) {
    float t    = closest_t(p0, p1, p2, p3, x);
    vec3 cp    = bspline_pos(p0, p1, p2, p3, t);
    vec3 d     = x - cp;
    float dist = sqrtf(dot(d, d));
    if (dist < 1e-7f) {
        out_grad_mu = {0.f, 1.f, 0.f};
        return r;
    }
    out_grad_mu = d * (1.f / dist); // outward radial unit vector = ∇μ
    return r - dist;
}

// ---------------------------------------------------------------------------
// Sparse Convolution Noise — Gaussian (squared-exponential) kernel
// Kernel: h(r) = exp(−r²/(2l²)),  l = cs/3
// Truncated to zero for r > cs (cell support radius).
// Variance (3D, infinite-integral approximation, truncation error ~1%):
//   Var = (N/cs³) · π^(3/2) · l³  = N · π^(3/2) / 27
// ---------------------------------------------------------------------------

static constexpr int SC_IMPULSES = 10; // impulses per cell

// sqrt(π) ≈ 1.7724538509
static constexpr float SC_KERNEL_VAR = (float)SC_IMPULSES * 3.14159265f * 1.7724538509f / 27.f;

// Truncated Gaussian kernel: h(d²) = exp(−d²·inv_l2/2), inv_l2 = 1/l² = 9/cs²
// Returns 0 implicitly — caller must check d2 < cs2 before calling.
inline float sc_kernel(float d2, float inv_l2) {
    return expf(-d2 * inv_l2 * 0.5f);
}

// κ(r²) = exp(−r²/(4l²))  — stationary covariance normalized to κ(0)=1
inline float kappa(float r2, float l2) {
    return expf(-r2 / (4.f * l2));
}
// κ''(0) = −1/(2l²)
inline float kappa_dd0(float l2) {
    return -1.f / (2.f * l2);
}

// Deterministic hash for cell (ix, iy, iz) and global seed
static uint32_t cell_hash(int32_t ix, int32_t iy, int32_t iz, uint32_t seed) {
    // FNV-1a-style mixing
    uint32_t h = seed ^ 2166136261u;
    h ^= (uint32_t)ix * 2654435761u;
    h ^= h >> 16;
    h *= 0x45d9f3bu;
    h ^= (uint32_t)iy * 2246822519u;
    h ^= h >> 16;
    h *= 0x45d9f3bu;
    h ^= (uint32_t)iz * 3266489917u;
    h ^= h >> 16;
    h *= 0x45d9f3bu;
    return h;
}

// LCG step: produces a float in [0, 1)
static float lcg_next(uint32_t &s) {
    s = s * 1664525u + 1013904223u;
    return (s >> 8) * (1.f / (1u << 24));
}

// Evaluate sparse convolution noise — value only (fast path, used during marching)
static float scnoise_val(vec3 p, float cell_size, uint32_t seed) {
    const float inv_cs   = 1.f / cell_size;
    const float cs2      = cell_size * cell_size;
    const float inv_l2   = 9.f * inv_cs * inv_cs; // 1/l², l = cs/3
    const float inv_norm = 1.f / sqrtf(SC_KERNEL_VAR);

    const int gx = (int)floorf(p.x * inv_cs);
    const int gy = (int)floorf(p.y * inv_cs);
    const int gz = (int)floorf(p.z * inv_cs);

    float vsum = 0.f;

    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz) {
                const int cx = gx + dx, cy = gy + dy, cz = gz + dz;
                uint32_t cs = cell_hash(cx, cy, cz, seed);

                for (int k = 0; k < SC_IMPULSES; ++k) {
                    const float xi = cell_size * ((float)cx + lcg_next(cs));
                    const float yi = cell_size * ((float)cy + lcg_next(cs));
                    const float zi = cell_size * ((float)cz + lcg_next(cs));
                    const float wi = (lcg_next(cs) < 0.5f) ? 1.f : -1.f;

                    const float rx = p.x - xi, ry = p.y - yi, rz = p.z - zi;
                    const float d2 = rx * rx + ry * ry + rz * rz;
                    if (d2 >= cs2)
                        continue;
                    vsum += wi * sc_kernel(d2, inv_l2);
                }
            }
    return vsum * inv_norm;
}

// Evaluate sparse convolution noise — value + gradient (used at hit point for normal)
static void
scnoise_val_grad(vec3 p, float cell_size, uint32_t seed, float &out_val, vec3 &out_grad) {
    const float inv_cs   = 1.f / cell_size;
    const float cs2      = cell_size * cell_size;
    const float inv_l2   = 9.f * inv_cs * inv_cs; // 1/l², l = cs/3
    const float inv_norm = 1.f / sqrtf(SC_KERNEL_VAR);

    const int gx = (int)floorf(p.x * inv_cs);
    const int gy = (int)floorf(p.y * inv_cs);
    const int gz = (int)floorf(p.z * inv_cs);

    float vsum = 0.f;
    vec3 gsum  = {0.f, 0.f, 0.f};

    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz) {
                const int cx = gx + dx, cy = gy + dy, cz = gz + dz;
                uint32_t cs = cell_hash(cx, cy, cz, seed);

                for (int k = 0; k < SC_IMPULSES; ++k) {
                    const float xi = cell_size * ((float)cx + lcg_next(cs));
                    const float yi = cell_size * ((float)cy + lcg_next(cs));
                    const float zi = cell_size * ((float)cz + lcg_next(cs));
                    const float wi = (lcg_next(cs) < 0.5f) ? 1.f : -1.f;

                    const float rx = p.x - xi, ry = p.y - yi, rz = p.z - zi;
                    const float d2 = rx * rx + ry * ry + rz * rz;
                    if (d2 >= cs2)
                        continue;

                    float kv  = sc_kernel(d2, inv_l2);
                    float dkc = -inv_l2 * kv; // ∂h/∂r_j = -(r_j/l²)·h → factor per-component
                    vsum += wi * kv;
                    gsum.x += wi * dkc * rx;
                    gsum.y += wi * dkc * ry;
                    gsum.z += wi * dkc * rz;
                }
            }
    out_val  = vsum * inv_norm;
    out_grad = gsum * inv_norm;
}

// ---------------------------------------------------------------------------
// Ray–AABB slab test
// ---------------------------------------------------------------------------
static bool ray_aabb(float ox,
                     float oy,
                     float oz,
                     float idx,
                     float idy,
                     float idz,
                     float tnear,
                     float tfar,
                     const RTCBounds &b,
                     float &t0,
                     float &t1) {
    float tx0 = (b.lower_x - ox) * idx, tx1 = (b.upper_x - ox) * idx;
    float ty0 = (b.lower_y - oy) * idy, ty1 = (b.upper_y - oy) * idy;
    float tz0 = (b.lower_z - oz) * idz, tz1 = (b.upper_z - oz) * idz;
    if (tx0 > tx1)
        std::swap(tx0, tx1);
    if (ty0 > ty1)
        std::swap(ty0, ty1);
    if (tz0 > tz1)
        std::swap(tz0, tz1);
    t0 = std::max({tx0, ty0, tz0, tnear});
    t1 = std::min({tx1, ty1, tz1, tfar});
    return t0 <= t1;
}

// ---------------------------------------------------------------------------
// AABB for one segment.
// Expanded by the noise amplitude so the ray always starts outside the
// noisy surface, guaranteeing f(t0) < 0 at AABB entry.
// ---------------------------------------------------------------------------
static RTCBounds segment_aabb(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float r, float amplitude) {
    // Expand by r + noise headroom so f < 0 at AABB boundary
    const float rr = r * (1.f + 3.f * amplitude);
    RTCBounds b;
    b.lower_x = std::min({p0.x, p1.x, p2.x, p3.x}) - rr;
    b.lower_y = std::min({p0.y, p1.y, p2.y, p3.y}) - rr;
    b.lower_z = std::min({p0.z, p1.z, p2.z, p3.z}) - rr;
    b.upper_x = std::max({p0.x, p1.x, p2.x, p3.x}) + rr;
    b.upper_y = std::max({p0.y, p1.y, p2.y, p3.y}) + rr;
    b.upper_z = std::max({p0.z, p1.z, p2.z, p3.z}) + rr;
    return b;
}

// ---------------------------------------------------------------------------
// Conditioning state — set by renderer before each rtcIntersect / rtcOccluded
// ---------------------------------------------------------------------------
thread_local CondState tl_cond = {{0.f, 0.f, 0.f}, 0.f, false};

// ---------------------------------------------------------------------------
// Shared march state for intersect and occluded
// ---------------------------------------------------------------------------
struct MarchResult {
    float t_hit;
    float ct;        // B-spline parameter (→ RTCHit.u, unchanged)
    vec3 normal;     // normalized ∇f at hit
    float u_val;     // pathwise value weight   (→ RTCHit.v)
    vec3 u_grad;     // pathwise gradient weights (→ tl_cond)
    float gz;        // ∂ψ/∂z at hit (for NEE)
    float noise_val; // ψ(x_s) raw
    bool hit;
};

static MarchResult march_segment(float ox,
                                 float oy,
                                 float oz,
                                 float dx,
                                 float dy,
                                 float dz,
                                 float tff,
                                 vec3 p0,
                                 vec3 p1,
                                 vec3 p2,
                                 vec3 p3,
                                 float r,
                                 const RTCBounds &aabb,
                                 float tnn,
                                 float amplitude, // = hair->amplitude * r  (world-space amplitude)
                                 uint32_t seed) {
    MarchResult result{};
    result.hit = false;

    float idx = (fabsf(dx) > 1e-9f) ? 1.f / dx : 1e30f * (dx >= 0.f ? 1.f : -1.f);
    float idy = (fabsf(dy) > 1e-9f) ? 1.f / dy : 1e30f * (dy >= 0.f ? 1.f : -1.f);
    float idz = (fabsf(dz) > 1e-9f) ? 1.f / dz : 1e30f * (dz >= 0.f ? 1.f : -1.f);

    float t0, t1;
    if (!ray_aabb(ox, oy, oz, idx, idy, idz, tnn, tff, aabb, t0, t1))
        return result;

    const float cell_size = r;                           // noise grid cell size = hair radius
    const float l2        = cell_size * cell_size / 9.f; // l = cs/3

    // Lipschitz-bounded step size (STEP 3)
    // L_ψ = sqrt(N_overlap) * (1/l) * exp(-0.5) / sqrt(Var)
    // N_overlap = SC_IMPULSES * 27 (worst-case overlapping cells)
    const float inv_norm    = 1.f / sqrtf(SC_KERNEL_VAR);
    const float inv_l       = 3.f / cell_size;
    const float N_overlap_f = (float)(SC_IMPULSES * 27);
    const float L_psi       = sqrtf(N_overlap_f) * inv_l * expf(-0.5f) * inv_norm;
    const float L_f         = 1.f + amplitude * L_psi;
    const float step_safe   = 1.f / L_f;

    const float span = t1 - t0;
    float step       = std::min(span / 4.f, step_safe);
    int nsteps       = (int)ceilf(span / step) + 1;
    nsteps           = std::min(nsteps, 256);

    // Value-only evaluation for marching — includes conditioning correction
    auto eval_f = [&](float t) -> float {
        vec3 x   = {ox + dx * t, oy + dy * t, oz + dz * t};
        float mu = mean_field_val(p0, p1, p2, p3, r, x);
        float f  = mu + amplitude * scnoise_val(x, cell_size, seed);

        if (tl_cond.valid) {
            // x0 = ray origin = (ox, oy, oz) = previous bounce hit point
            vec3 delta = {x.x - ox, x.y - oy, x.z - oz};
            float r2   = dot(delta, delta);
            float kv   = kappa(r2, l2);
            // ∇₂κ(x, x0) = −(x−x0)/(2l²) · κ(x, x0)  (gradient wrt x0)
            // dot(∇₂κ, u_grad) = −(1/(2l²)) · kv · dot(delta, u_grad)
            float grad_splat = (-1.f / (2.f * l2)) * kv * dot(delta, tl_cond.u_grad);
            f += kv * tl_cond.u_val + grad_splat;
        }

        return f;
    };

    float f_prev = eval_f(t0);
    int sign0    = (f_prev >= 0.f) ? 1 : -1;

    for (int s = 1; s <= nsteps; ++s) {
        float t_cur = t0 + s * step;
        float f_cur = eval_f(t_cur);
        int sign1   = (f_cur >= 0.f) ? 1 : -1;

        if (sign0 < 0 && sign1 > 0) {
            // Bisect between t_cur - step and t_cur
            float ta = t_cur - step;
            float tb = t_cur;
            for (int iter = 0; iter < 8; ++iter) {
                float tm = 0.5f * (ta + tb);
                float fm = eval_f(tm);
                if ((fm >= 0.f ? 1 : -1) == sign0)
                    ta = tm;
                else
                    tb = tm;
            }
            float t_hit = 0.5f * (ta + tb);
            if (t_hit >= tff)
                return result;

            // Normal = ∇f at hit = ∇μ + amplitude · ∇noise (full gradient only here)
            vec3 x_hit = {ox + dx * t_hit, oy + dy * t_hit, oz + dz * t_hit};
            vec3 grad_mu;
            mean_field(p0, p1, p2, p3, r, x_hit, grad_mu);
            float noise_val;
            vec3 noise_grad;
            scnoise_val_grad(x_hit, cell_size, seed, noise_val, noise_grad);

            vec3 ng      = grad_mu + noise_grad * amplitude;
            float ng_len = sqrtf(dot(ng, ng));
            if (ng_len < 1e-7f)
                ng = grad_mu;
            else
                ng = ng * (1.f / ng_len);

            // Curve parameter for tangent reconstruction
            float ct = closest_t(p0, p1, p2, p3, x_hit);

            // Pathwise conditioning weights (STEP 5.3)
            float u_val = -noise_val;
            vec3 u_grad = noise_grad * (-2.f * l2);

            // Directional noise derivative along ray direction (for NEE)
            float gz = dot(noise_grad, {dx, dy, dz});

            result.hit       = true;
            result.t_hit     = t_hit;
            result.ct        = ct;
            result.normal    = ng;
            result.u_val     = u_val;
            result.u_grad    = u_grad;
            result.gz        = gz;
            result.noise_val = noise_val;
            return result;
        }

        f_prev = f_cur;
    }

    return result;
}

// ---------------------------------------------------------------------------
// Embree callbacks
// ---------------------------------------------------------------------------

static void bounds_cb(const RTCBoundsFunctionArguments *args) {
    const HairData *hair = static_cast<const HairData *>(args->geometryUserPtr);
    int base             = hair->indices[args->primID];
    float r              = std::max({hair->vertices[base].w,
                                     hair->vertices[base + 1].w,
                                     hair->vertices[base + 2].w,
                                     hair->vertices[base + 3].w});
    vec3 p0              = from_f4(hair->vertices[base + 0]);
    vec3 p1              = from_f4(hair->vertices[base + 1]);
    vec3 p2              = from_f4(hair->vertices[base + 2]);
    vec3 p3              = from_f4(hair->vertices[base + 3]);
    *args->bounds_o      = segment_aabb(p0, p1, p2, p3, r, hair->amplitude);
}

static void intersect_cb(const RTCIntersectFunctionNArguments *args) {
    const HairData *hair    = static_cast<const HairData *>(args->geometryUserPtr);
    const unsigned int N    = args->N;
    const unsigned int prim = args->primID;
    const unsigned int geom = args->geomID;

    RTCRayN *rayN = RTCRayHitN_RayN(args->rayhit, N);
    RTCHitN *hitN = RTCRayHitN_HitN(args->rayhit, N);

    int base = hair->indices[prim];
    float r  = std::max({hair->vertices[base].w,
                         hair->vertices[base + 1].w,
                         hair->vertices[base + 2].w,
                         hair->vertices[base + 3].w});
    vec3 p0  = from_f4(hair->vertices[base + 0]);
    vec3 p1  = from_f4(hair->vertices[base + 1]);
    vec3 p2  = from_f4(hair->vertices[base + 2]);
    vec3 p3  = from_f4(hair->vertices[base + 3]);

    RTCBounds aabb        = segment_aabb(p0, p1, p2, p3, r, hair->amplitude);
    const float amplitude = hair->amplitude * r;

    for (unsigned int i = 0; i < N; ++i) {
        if (!args->valid[i])
            continue;

        float ox  = RTCRayN_org_x(rayN, N, i);
        float oy  = RTCRayN_org_y(rayN, N, i);
        float oz  = RTCRayN_org_z(rayN, N, i);
        float ddx = RTCRayN_dir_x(rayN, N, i);
        float ddy = RTCRayN_dir_y(rayN, N, i);
        float ddz = RTCRayN_dir_z(rayN, N, i);
        float tnn = RTCRayN_tnear(rayN, N, i);
        float tff = RTCRayN_tfar(rayN, N, i);

        MarchResult res = march_segment(
            ox, oy, oz, ddx, ddy, ddz, tff, p0, p1, p2, p3, r, aabb, tnn, amplitude, hair->seed);
        if (!res.hit)
            continue;

        RTCRayN_tfar(rayN, N, i)   = res.t_hit;
        RTCHitN_Ng_x(hitN, N, i)   = res.normal.x;
        RTCHitN_Ng_y(hitN, N, i)   = res.normal.y;
        RTCHitN_Ng_z(hitN, N, i)   = res.normal.z;
        RTCHitN_u(hitN, N, i)      = res.ct;    // B-spline parameter (unchanged)
        RTCHitN_v(hitN, N, i)      = res.u_val; // pathwise value weight
        RTCHitN_primID(hitN, N, i) = prim;
        RTCHitN_geomID(hitN, N, i) = geom;
        tl_cond = {res.u_grad, res.u_val, true}; // gradient weights for next bounce
    }
}

// ---------------------------------------------------------------------------

static void occluded_cb(const RTCOccludedFunctionNArguments *args) {
    const HairData *hair    = static_cast<const HairData *>(args->geometryUserPtr);
    const unsigned int N    = args->N;
    const unsigned int prim = args->primID;

    RTCRayN *rayN = args->ray;

    int base = hair->indices[prim];
    float r  = std::max({hair->vertices[base].w,
                         hair->vertices[base + 1].w,
                         hair->vertices[base + 2].w,
                         hair->vertices[base + 3].w});
    vec3 p0  = from_f4(hair->vertices[base + 0]);
    vec3 p1  = from_f4(hair->vertices[base + 1]);
    vec3 p2  = from_f4(hair->vertices[base + 2]);
    vec3 p3  = from_f4(hair->vertices[base + 3]);

    RTCBounds aabb        = segment_aabb(p0, p1, p2, p3, r, hair->amplitude);
    const float amplitude = hair->amplitude * r;

    for (unsigned int i = 0; i < N; ++i) {
        if (!args->valid[i])
            continue;

        float ox  = RTCRayN_org_x(rayN, N, i);
        float oy  = RTCRayN_org_y(rayN, N, i);
        float oz  = RTCRayN_org_z(rayN, N, i);
        float ddx = RTCRayN_dir_x(rayN, N, i);
        float ddy = RTCRayN_dir_y(rayN, N, i);
        float ddz = RTCRayN_dir_z(rayN, N, i);
        float tnn = RTCRayN_tnear(rayN, N, i);
        float tff = RTCRayN_tfar(rayN, N, i);

        MarchResult res = march_segment(
            ox, oy, oz, ddx, ddy, ddz, tff, p0, p1, p2, p3, r, aabb, tnn, amplitude, hair->seed);
        if (res.hit) {
            // Signal occlusion: Embree checks tfar < 0
            // Use bit-cast to avoid UB under -ffast-math
            static const uint32_t neg_inf_bits = 0xFF800000u;
            float neg_inf;
            __builtin_memcpy(&neg_inf, &neg_inf_bits, 4);
            RTCRayN_tfar(rayN, N, i) = neg_inf;
        }
    }
}

// ---------------------------------------------------------------------------

GpisHitInfo extract_hit_info(const RTCRayHit &rh, const HairData &hair) {
    const float r         = hair.vertices[0].w; // approximate — caller can refine
    const float cell_size = r;
    const float l2        = cell_size * cell_size / 9.f;

    vec3 p = {rh.ray.org_x + rh.ray.dir_x * rh.ray.tfar,
              rh.ray.org_y + rh.ray.dir_y * rh.ray.tfar,
              rh.ray.org_z + rh.ray.dir_z * rh.ray.tfar};

    // Recover noise gradient from hit normal and μ gradient
    // ng_normalized = normalize(grad_mu + amplitude * noise_grad)
    // We store Ng un-normalized (just the direction); compute grad_mu at p.
    // For NEE, caller needs the raw noise gradient — reconstruct via:
    //   noise_grad ≈ (ng - grad_mu) / amplitude  (if amplitude > 0)
    // Use the pre-computed gz (∂ψ/∂z) stored via tl_cond or caller-computed.

    GpisHitInfo info;
    info.p    = p;
    info.l2   = l2;
    info.kdd0 = kappa_dd0(l2);
    info.gz   = 0.f; // caller should fill from saved MarchResult.gz

    // ∇μ at hit point
    float ct = rh.hit.u; // B-spline parameter
    (void)ct;
    // Ng is the un-normalized surface normal stored in the hit record
    info.grad_mu = normalize({rh.hit.Ng_x, rh.hit.Ng_y, rh.hit.Ng_z});

    return info;
}

// ---------------------------------------------------------------------------

unsigned add_user_hair(RTCDevice device, RTCScene scene, const HairData &hair) {
    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
    rtcSetGeometryUserPrimitiveCount(geom, (unsigned)hair.indices.size());
    rtcSetGeometryUserData(geom, const_cast<HairData *>(&hair));
    rtcSetGeometryBoundsFunction(geom, bounds_cb, nullptr);
    rtcSetGeometryIntersectFunction(geom, intersect_cb);
    rtcSetGeometryOccludedFunction(geom, occluded_cb);
    rtcCommitGeometry(geom);
    unsigned id = rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    return id;
}

} // namespace m3hair
