#include "gpis_geo.h"

#include "gpis_nee.h"
#include "math.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace m3hair {

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

static constexpr int SC_IMPULSES = 10;
static constexpr float SC_KERNEL_VAR = (float)SC_IMPULSES * 3.14159265f * 1.7724538509f / 27.f;

inline float sc_kernel(float d2, float inv_l2) {
    return expf(-d2 * inv_l2 * 0.5f);
}

inline float kappa(float r2, float l2) {
    return expf(-r2 / (4.f * l2));
}

inline float kappa_dd0(float l2) {
    return -1.f / (2.f * l2);
}

static uint32_t cell_hash(int32_t ix, int32_t iy, int32_t iz, uint32_t seed) {
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

static float lcg_next(uint32_t &s) {
    s = s * 1664525u + 1013904223u;
    return (s >> 8) * (1.f / (1u << 24));
}

static float scnoise_val(vec3 p, float cell_size, uint32_t seed) {
    const float inv_cs   = 1.f / cell_size;
    const float cs2      = cell_size * cell_size;
    const float inv_l2   = 9.f * inv_cs * inv_cs;
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

static void
scnoise_val_grad(vec3 p, float cell_size, uint32_t seed, float &out_val, vec3 &out_grad) {
    const float inv_cs   = 1.f / cell_size;
    const float cs2      = cell_size * cell_size;
    const float inv_l2   = 9.f * inv_cs * inv_cs;
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
                    float dkc = -inv_l2 * kv;
                    vsum += wi * kv;
                    gsum.x += wi * dkc * rx;
                    gsum.y += wi * dkc * ry;
                    gsum.z += wi * dkc * rz;
                }
            }
    out_val  = vsum * inv_norm;
    out_grad = gsum * inv_norm;
}

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

static RTCBounds segment_aabb(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float r, float amplitude) {
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

thread_local CondState tl_cond = {{0.f, 0.f, 0.f}, 0.f, false};

struct MarchResult {
    float t_hit;
    float ct;
    vec3 normal;
    float u_val;
    vec3 u_grad;
    float gz;
    float noise_val;
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
                                 float amplitude,
                                 float cell_size,
                                 uint32_t seed) {
    MarchResult result{};

    float idx = (fabsf(dx) > 1e-9f) ? 1.f / dx : 1e30f * (dx >= 0.f ? 1.f : -1.f);
    float idy = (fabsf(dy) > 1e-9f) ? 1.f / dy : 1e30f * (dy >= 0.f ? 1.f : -1.f);
    float idz = (fabsf(dz) > 1e-9f) ? 1.f / dz : 1e30f * (dz >= 0.f ? 1.f : -1.f);

    float t0, t1;
    if (!ray_aabb(ox, oy, oz, idx, idy, idz, tnn, tff, aabb, t0, t1))
        return result;

    const vec3 O = {ox, oy, oz};
    const vec3 D = {dx, dy, dz};

    constexpr int N = 16;
    vec3 caps[N + 1];
    for (int i = 0; i <= N; ++i)
        caps[i] = bspline_pos(p0, p1, p2, p3, i / (float)N);

    const float L_noise = (amplitude > 0.f && cell_size > 1e-8f)
                        ? SC_IMPULSES * 27.f * 0.60653f * (3.f / cell_size) / sqrtf(SC_KERNEL_VAR)
                        : 0.f;
    const float L_total = 1.f + amplitude * L_noise;

    const float hit_eps  = r * 1e-3f;
    const float min_step = r * 1e-5f;
    constexpr int MAX_STEPS = 256;

    float t = t0;

    for (int step = 0; step < MAX_STEPS; ++step) {
        if (t > t1) break;

        const vec3 Q = O + D * t;

        float min_dist_sq = 1e30f;
        vec3  closest_v   = {};

        for (int i = 0; i < N; ++i) {
            const vec3  A    = caps[i];
            const vec3  AB   = caps[i + 1] - A;
            const float len2 = dot(AB, AB);
            const float s    = (len2 > 1e-12f)
                             ? std::clamp(dot(Q - A, AB) / len2, 0.f, 1.f)
                             : 0.f;
            const vec3  v    = Q - (A + AB * s);
            const float d2   = dot(v, v);
            if (d2 < min_dist_sq) { min_dist_sq = d2; closest_v = v; }
        }

        const float min_dist = sqrtf(min_dist_sq);

        float noise_val  = 0.f;
        vec3  noise_grad = {};
        if (amplitude > 0.f)
            scnoise_val_grad(Q, cell_size, seed, noise_val, noise_grad);

        const float phi = min_dist - r - amplitude * noise_val;

        if (phi <= hit_eps) {
            const vec3 n_cap  = (min_dist > 1e-8f)
                              ? closest_v * (1.f / min_dist)
                              : vec3{0.f, 1.f, 0.f};
            const vec3 normal = normalize(n_cap - noise_grad * amplitude);
            return {t, 0.f, normal, 0.f, {0.f, 0.f, 0.f}, 0.f, 0.f, true};
        }

        t += std::max(phi / L_total, min_step);
    }

    return result;
}

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
            ox, oy, oz, ddx, ddy, ddz, tff, p0, p1, p2, p3, r, aabb, tnn,
            amplitude, hair->cell_size, hair->seed);
        if (!res.hit)
            continue;

        RTCRayN_tfar(rayN, N, i)   = res.t_hit;
        RTCHitN_Ng_x(hitN, N, i)   = res.normal.x;
        RTCHitN_Ng_y(hitN, N, i)   = res.normal.y;
        RTCHitN_Ng_z(hitN, N, i)   = res.normal.z;
        RTCHitN_u(hitN, N, i)      = res.ct;
        RTCHitN_v(hitN, N, i)      = res.u_val;
        RTCHitN_primID(hitN, N, i) = prim;
        RTCHitN_geomID(hitN, N, i) = geom;
        tl_cond = {res.u_grad, res.u_val, true};
    }
}

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
            ox, oy, oz, ddx, ddy, ddz, tff, p0, p1, p2, p3, r, aabb, tnn,
            amplitude, hair->cell_size, hair->seed);
        if (res.hit) {
            // bit-cast to avoid UB under -ffast-math
            static const uint32_t neg_inf_bits = 0xFF800000u;
            float neg_inf;
            __builtin_memcpy(&neg_inf, &neg_inf_bits, 4);
            RTCRayN_tfar(rayN, N, i) = neg_inf;
        }
    }
}

GpisHitInfo extract_hit_info(const RTCRayHit &rh, const HairData &hair) {
    const float r  = hair.vertices[0].w;
    const float l2 = r * r / 9.f;

    return {{rh.ray.org_x + rh.ray.dir_x * rh.ray.tfar,
             rh.ray.org_y + rh.ray.dir_y * rh.ray.tfar,
             rh.ray.org_z + rh.ray.dir_z * rh.ray.tfar},
            normalize({rh.hit.Ng_x, rh.hit.Ng_y, rh.hit.Ng_z}),
            0.f,
            l2,
            kappa_dd0(l2)};
}

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
