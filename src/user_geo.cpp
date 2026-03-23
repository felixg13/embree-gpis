// GPIS (Gaussian Process Implicit Surface) user geometry for Embree 4.
//
// Each B-spline segment is one user primitive. The intersect callback:
//   1. Marches the ray through the segment AABB (step ≈ radius)
//   2. Samples n=8 GP points along the ray
//   3. Builds covariance K with squared-exponential kernel (l = radius)
//   4. Draws a realization f = μ + chol(K)*z  (z seeded from primID)
//   5. Bisects the first sign change to find t_hit
//   6. Normal = outward radial from closest curve point (gradient of μ)
//
// The z vector is fixed per primID → same surface shape regardless of ray,
// giving temporal coherence with no extra storage.

#include "user_geo.h"
#include "math.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace m3hair {

static constexpr int NGP = 8;   // GP sample points per ray

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

static vec3 from_f4(const float4& v) { return {v.x, v.y, v.z}; }

// Cubic uniform B-spline position at t ∈ [0,1]
static vec3 bspline_pos(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float t) {
    float t2 = t*t, t3 = t2*t;
    float b0 = (-t3 + 3*t2 - 3*t + 1) * (1.f/6.f);
    float b1 = ( 3*t3 - 6*t2      + 4) * (1.f/6.f);
    float b2 = (-3*t3 + 3*t2 + 3*t+ 1) * (1.f/6.f);
    float b3 = t3 * (1.f/6.f);
    return p0*b0 + p1*b1 + p2*b2 + p3*b3;
}

// Approximate closest t on B-spline to point x via 12-sample search
static float closest_t(vec3 p0, vec3 p1, vec3 p2, vec3 p3, vec3 x) {
    float best_t = 0.f, best_d2 = 1e30f;
    for (int i = 0; i <= 12; ++i) {
        float t  = i / 12.f;
        vec3  pt = bspline_pos(p0, p1, p2, p3, t);
        float d2 = dot(pt - x, pt - x);
        if (d2 < best_d2) { best_d2 = d2; best_t = t; }
    }
    // One Newton refinement step (approx)
    float eps = 1.f/24.f;
    for (int iter = 0; iter < 2; ++iter) {
        vec3 pa = bspline_pos(p0, p1, p2, p3, std::max(0.f, best_t - eps));
        vec3 pb = bspline_pos(p0, p1, p2, p3, std::min(1.f, best_t + eps));
        float da = dot(pa - x, pa - x), db = dot(pb - x, pb - x);
        if (da < best_d2) { best_d2 = da; best_t = std::max(0.f, best_t-eps); }
        if (db < best_d2) { best_d2 = db; best_t = std::min(1.f, best_t+eps); }
        eps *= 0.5f;
    }
    return best_t;
}

// Mean function: r - dist_to_bspline  (positive inside cylinder, negative outside)
static float mean_field(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float r, vec3 x) {
    float t  = closest_t(p0, p1, p2, p3, x);
    vec3  cp = bspline_pos(p0, p1, p2, p3, t);
    return r - sqrtf(dot(x - cp, x - cp));
}

// Squared-exponential kernel
static float kernel(vec3 a, vec3 b, float inv2l2) {
    float d2 = dot(a - b, a - b);
    return expf(-d2 * inv2l2);
}

// ---------------------------------------------------------------------------
// Inline 8×8 Cholesky (lower triangular, in-place)
// ---------------------------------------------------------------------------
static void cholesky8(float L[NGP][NGP]) {
    for (int j = 0; j < NGP; ++j) {
        float s = L[j][j];
        for (int k = 0; k < j; ++k) s -= L[j][k]*L[j][k];
        L[j][j] = (s > 0.f) ? sqrtf(s) : 0.f;
        float inv = (L[j][j] > 1e-9f) ? 1.f/L[j][j] : 0.f;
        for (int i = j+1; i < NGP; ++i) {
            float r = L[i][j];
            for (int k = 0; k < j; ++k) r -= L[i][k]*L[j][k];
            L[i][j] = r * inv;
        }
        for (int i = 0; i < j; ++i) L[j][i] = 0.f;   // zero upper triangle
    }
}

// ---------------------------------------------------------------------------
// Deterministic Gaussian z[8] from primID seed (Box-Muller)
// ---------------------------------------------------------------------------
static uint32_t xorshift(uint32_t s) {
    s ^= s << 13; s ^= s >> 17; s ^= s << 5; return s;
}
static float u2f(uint32_t s) { return (s >> 8) * (1.f/(1u<<24)); }

static void gen_z(uint32_t seed, float z[NGP]) {
    seed = xorshift(seed * 2654435761u + 1u);
    for (int i = 0; i < NGP/2; ++i) {
        seed = xorshift(seed); float u1 = std::max(u2f(seed), 1e-7f);
        seed = xorshift(seed); float u2 = u2f(seed);
        float r = sqrtf(-2.f * logf(u1));
        float th = 6.28318530f * u2;
        z[2*i]   = r * cosf(th);
        z[2*i+1] = r * sinf(th);
    }
}

// ---------------------------------------------------------------------------
// Ray–AABB slab test; returns false if no overlap with [tnear,tfar]
// ---------------------------------------------------------------------------
static bool ray_aabb(float ox, float oy, float oz,
                     float idx, float idy, float idz,
                     float tnear, float tfar,
                     const RTCBounds& b,
                     float& t0, float& t1)
{
    float tx0 = (b.lower_x - ox) * idx, tx1 = (b.upper_x - ox) * idx;
    float ty0 = (b.lower_y - oy) * idy, ty1 = (b.upper_y - oy) * idy;
    float tz0 = (b.lower_z - oz) * idz, tz1 = (b.upper_z - oz) * idz;
    if (tx0 > tx1) std::swap(tx0, tx1);
    if (ty0 > ty1) std::swap(ty0, ty1);
    if (tz0 > tz1) std::swap(tz0, tz1);
    t0 = std::max({tx0, ty0, tz0, tnear});
    t1 = std::min({tx1, ty1, tz1, tfar});
    return t0 <= t1;
}

// ---------------------------------------------------------------------------
// AABB for a single segment (convex hull of 4 CPs ± max_radius)
// ---------------------------------------------------------------------------
static RTCBounds segment_aabb(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float r) {
    RTCBounds b;
    b.lower_x = std::min({p0.x,p1.x,p2.x,p3.x}) - r;
    b.lower_y = std::min({p0.y,p1.y,p2.y,p3.y}) - r;
    b.lower_z = std::min({p0.z,p1.z,p2.z,p3.z}) - r;
    b.upper_x = std::max({p0.x,p1.x,p2.x,p3.x}) + r;
    b.upper_y = std::max({p0.y,p1.y,p2.y,p3.y}) + r;
    b.upper_z = std::max({p0.z,p1.z,p2.z,p3.z}) + r;
    return b;
}

// ---------------------------------------------------------------------------
// Embree callbacks
// ---------------------------------------------------------------------------

static void bounds_cb(const RTCBoundsFunctionArguments* args) {
    const HairData* hair = static_cast<const HairData*>(args->geometryUserPtr);
    int base = hair->indices[args->primID];
    float r = std::max({hair->vertices[base].w, hair->vertices[base+1].w,
                        hair->vertices[base+2].w, hair->vertices[base+3].w});
    vec3 p0 = from_f4(hair->vertices[base+0]);
    vec3 p1 = from_f4(hair->vertices[base+1]);
    vec3 p2 = from_f4(hair->vertices[base+2]);
    vec3 p3 = from_f4(hair->vertices[base+3]);
    *args->bounds_o = segment_aabb(p0, p1, p2, p3, r);
}

static void intersect_cb(const RTCIntersectFunctionNArguments* args) {
    const HairData* hair = static_cast<const HairData*>(args->geometryUserPtr);
    const unsigned int N    = args->N;
    const unsigned int prim = args->primID;
    const unsigned int geom = args->geomID;

    RTCRayN* rayN = RTCRayHitN_RayN(args->rayhit, N);
    RTCHitN* hitN = RTCRayHitN_HitN(args->rayhit, N);

    int   base = hair->indices[prim];
    float r    = std::max({hair->vertices[base].w, hair->vertices[base+1].w,
                           hair->vertices[base+2].w, hair->vertices[base+3].w});
    vec3 p0 = from_f4(hair->vertices[base+0]);
    vec3 p1 = from_f4(hair->vertices[base+1]);
    vec3 p2 = from_f4(hair->vertices[base+2]);
    vec3 p3 = from_f4(hair->vertices[base+3]);

    // Build AABB once — same for all rays in the packet
    RTCBounds aabb = segment_aabb(p0, p1, p2, p3, r);

    // Fixed z per primitive (deterministic noise = fixed surface shape)
    float z[NGP];
    gen_z(prim, z);

    // Kernel length scale and amplitude
    const float sigma2  = (r*r) * 0.04f;   // amplitude² = (0.2*r)²
    const float inv2l2  = 1.f / (2.f * r*r);
    const float nugget  = sigma2 * 1e-3f;   // numerical stabiliser

    for (unsigned int i = 0; i < N; ++i) {
        if (!args->valid[i]) continue;

        float ox  = RTCRayN_org_x(rayN, N, i);
        float oy  = RTCRayN_org_y(rayN, N, i);
        float oz  = RTCRayN_org_z(rayN, N, i);
        float dx  = RTCRayN_dir_x(rayN, N, i);
        float dy  = RTCRayN_dir_y(rayN, N, i);
        float dz  = RTCRayN_dir_z(rayN, N, i);
        float tnn = RTCRayN_tnear(rayN, N, i);
        float tff = RTCRayN_tfar (rayN, N, i);

        float idx = (fabsf(dx) > 1e-9f) ? 1.f/dx : 1e30f * (dx >= 0.f ? 1.f : -1.f);
        float idy = (fabsf(dy) > 1e-9f) ? 1.f/dy : 1e30f * (dy >= 0.f ? 1.f : -1.f);
        float idz = (fabsf(dz) > 1e-9f) ? 1.f/dz : 1e30f * (dz >= 0.f ? 1.f : -1.f);

        float t0, t1;
        if (!ray_aabb(ox,oy,oz, idx,idy,idz, tnn,tff, aabb, t0,t1)) continue;

        // --- sample NGP points uniformly inside the AABB segment ---
        float step = (t1 - t0) / (NGP - 1);
        vec3  pts[NGP];
        float mu[NGP], f[NGP];

        for (int j = 0; j < NGP; ++j) {
            float t = t0 + j * step;
            pts[j] = { ox + dx*t, oy + dy*t, oz + dz*t };
            mu[j]  = mean_field(p0, p1, p2, p3, r, pts[j]);
        }

        // --- covariance matrix K (symmetric, NGP×NGP) ---
        float K[NGP][NGP];
        for (int a = 0; a < NGP; ++a)
            for (int b = a; b < NGP; ++b) {
                float k = sigma2 * kernel(pts[a], pts[b], inv2l2);
                if (a == b) k += nugget;
                K[a][b] = K[b][a] = k;
            }

        // --- Cholesky and realise f = mu + L*z ---
        cholesky8(K);   // K is overwritten with L
        for (int a = 0; a < NGP; ++a) {
            float eps = 0.f;
            for (int b = 0; b <= a; ++b) eps += K[a][b] * z[b];
            f[a] = mu[a] + eps;
        }

        // --- find first sign change ---
        int sc = -1;
        for (int j = 0; j < NGP-1; ++j) {
            if (f[j] * f[j+1] < 0.f) { sc = j; break; }
        }
        if (sc < 0) continue;   // no intersection with this segment

        // --- bisect to refine t_hit ---
        float ta = t0 + sc * step, fa = f[sc];
        float tb = t0 + (sc+1) * step;
        for (int iter = 0; iter < 8; ++iter) {
            float tm = 0.5f*(ta+tb);
            vec3  pm = { ox+dx*tm, oy+dy*tm, oz+dz*tm };
            float fm = mean_field(p0, p1, p2, p3, r, pm);
            // approximate GP epsilon at midpoint (interpolate linearly)
            float alpha = (tm - ta) / (tb - ta);
            fm += (f[sc] * (1.f-alpha) + f[sc+1] * alpha) - (mu[sc]*(1.f-alpha) + mu[sc+1]*alpha);
            if (fa * fm < 0.f) { tb = tm; }
            else               { ta = tm; fa = fm; }
        }
        float t_hit = 0.5f*(ta+tb);
        if (t_hit >= tff) continue;   // behind existing hit

        // --- normal: outward radial from closest curve point ---
        vec3  p_hit = { ox+dx*t_hit, oy+dy*t_hit, oz+dz*t_hit };
        float ct = closest_t(p0, p1, p2, p3, p_hit);
        vec3  cp = bspline_pos(p0, p1, p2, p3, ct);
        vec3  ng = p_hit - cp;
        float ng_len = sqrtf(dot(ng, ng));
        if (ng_len < 1e-7f) ng = {0.f, 1.f, 0.f};
        else                ng = ng / ng_len;

        // --- write hit record ---
        RTCRayN_tfar(rayN, N, i) = t_hit;
        RTCHitN_Ng_x(hitN, N, i) = ng.x;
        RTCHitN_Ng_y(hitN, N, i) = ng.y;
        RTCHitN_Ng_z(hitN, N, i) = ng.z;
        RTCHitN_u(hitN, N, i) = ct;          // curve parameter at hit
        RTCHitN_v(hitN, N, i) = 0.f;
        RTCHitN_primID(hitN, N, i) = prim;
        RTCHitN_geomID(hitN, N, i) = geom;
    }
}

// ---------------------------------------------------------------------------
// Occluded callback — same raymarching logic, just signals hit by setting tfar=-inf
// ---------------------------------------------------------------------------

static void occluded_cb(const RTCOccludedFunctionNArguments* args) {
    const HairData* hair = static_cast<const HairData*>(args->geometryUserPtr);
    const unsigned int N    = args->N;
    const unsigned int prim = args->primID;

    RTCRayN* rayN = args->ray;

    int   base = hair->indices[prim];
    float r    = std::max({hair->vertices[base].w, hair->vertices[base+1].w,
                           hair->vertices[base+2].w, hair->vertices[base+3].w});
    vec3 p0 = from_f4(hair->vertices[base+0]);
    vec3 p1 = from_f4(hair->vertices[base+1]);
    vec3 p2 = from_f4(hair->vertices[base+2]);
    vec3 p3 = from_f4(hair->vertices[base+3]);

    RTCBounds aabb = segment_aabb(p0, p1, p2, p3, r);

    float z[NGP];
    gen_z(prim, z);

    const float sigma2 = (r*r) * 0.04f;
    const float inv2l2 = 1.f / (2.f * r*r);
    const float nugget = sigma2 * 1e-3f;

    for (unsigned int i = 0; i < N; ++i) {
        if (!args->valid[i]) continue;

        float ox  = RTCRayN_org_x(rayN, N, i);
        float oy  = RTCRayN_org_y(rayN, N, i);
        float oz  = RTCRayN_org_z(rayN, N, i);
        float dx  = RTCRayN_dir_x(rayN, N, i);
        float dy  = RTCRayN_dir_y(rayN, N, i);
        float dz  = RTCRayN_dir_z(rayN, N, i);
        float tnn = RTCRayN_tnear(rayN, N, i);
        float tff = RTCRayN_tfar (rayN, N, i);

        float idx = (fabsf(dx) > 1e-9f) ? 1.f/dx : 1e30f * (dx >= 0.f ? 1.f : -1.f);
        float idy = (fabsf(dy) > 1e-9f) ? 1.f/dy : 1e30f * (dy >= 0.f ? 1.f : -1.f);
        float idz = (fabsf(dz) > 1e-9f) ? 1.f/dz : 1e30f * (dz >= 0.f ? 1.f : -1.f);

        float t0, t1;
        if (!ray_aabb(ox,oy,oz, idx,idy,idz, tnn,tff, aabb, t0,t1)) continue;

        float step = (t1 - t0) / (NGP - 1);
        vec3  pts[NGP];
        float mu[NGP], f[NGP];

        for (int j = 0; j < NGP; ++j) {
            float t = t0 + j * step;
            pts[j] = { ox + dx*t, oy + dy*t, oz + dz*t };
            mu[j]  = mean_field(p0, p1, p2, p3, r, pts[j]);
        }

        float K[NGP][NGP];
        for (int a = 0; a < NGP; ++a)
            for (int b = a; b < NGP; ++b) {
                float k = sigma2 * kernel(pts[a], pts[b], inv2l2);
                if (a == b) k += nugget;
                K[a][b] = K[b][a] = k;
            }

        cholesky8(K);
        for (int a = 0; a < NGP; ++a) {
            float eps = 0.f;
            for (int b = 0; b <= a; ++b) eps += K[a][b] * z[b];
            f[a] = mu[a] + eps;
        }

        int sc = -1;
        for (int j = 0; j < NGP-1; ++j) {
            if (f[j] * f[j+1] < 0.f) { sc = j; break; }
        }
        if (sc < 0) continue;

        // Signal occlusion
        RTCRayN_tfar(rayN, N, i) = -INFINITY;
    }
}

// ---------------------------------------------------------------------------

unsigned add_user_hair(RTCDevice device, RTCScene scene, const HairData& hair) {
    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
    rtcSetGeometryUserPrimitiveCount(geom, (unsigned)hair.indices.size());
    rtcSetGeometryUserData(geom, const_cast<HairData*>(&hair));
    rtcSetGeometryBoundsFunction(geom, bounds_cb, nullptr);
    rtcSetGeometryIntersectFunction(geom, intersect_cb);
    rtcSetGeometryOccludedFunction(geom, occluded_cb);
    rtcCommitGeometry(geom);
    unsigned id = rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    return id;
}

} // namespace m3hair
