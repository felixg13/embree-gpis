#include "gpis_geo.h"

#include "gpis_nee.h"
#include <cmath>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>

namespace m3hair {

static Vec3 from_f4(const Float4 &v) {
    return {v.x, v.y, v.z};
}

static Vec3 bspline_pos(Vec3 p0, Vec3 p1, Vec3 p2, Vec3 p3, float t) {
    float t2 = t * t;
    float t3 = t2 * t;
    float b0 = (-t3 + (3 * t2) - (3 * t) + 1) * (1.F / 6.F);
    float b1 = ((3 * t3) - (6 * t2) + 4) * (1.F / 6.F);
    float b2 = ((-3 * t3) + (3 * t2) + (3 * t) + 1) * (1.F / 6.F);
    float b3 = t3 * (1.F / 6.F);
    return p0 * b0 + p1 * b1 + p2 * b2 + p3 * b3;
}

static constexpr int SC_IMPULSES = 10;
static constexpr float SC_KERNEL_VAR = static_cast<float>(SC_IMPULSES) * std::numbers::pi_v<float> * 1.7724538509F / 27.F;

inline float sc_kernel(float d2, float inv_l2) {
    return expf(-d2 * inv_l2 * 0.5F);
}

inline float kappa(float r2, float l2) {
    return expf(-r2 / (4.F * l2));
}

inline float kappa_dd0(float l2) {
    return -1.F / (2.F * l2);
}

static uint32_t cell_hash(int32_t ix, int32_t iy, int32_t iz, uint32_t seed) {
    uint32_t h = seed ^ 2166136261U;
    h ^= static_cast<uint32_t>(ix) * 2654435761U;
    h ^= h >> 16;
    h *= 0x45d9f3bU;
    h ^= static_cast<uint32_t>(iy) * 2246822519U;
    h ^= h >> 16;
    h *= 0x45d9f3bU;
    h ^= static_cast<uint32_t>(iz) * 3266489917U;
    h ^= h >> 16;
    h *= 0x45d9f3bU;
    return h;
}

static float lcg_next(uint32_t &s) {
    s = (s * 1664525U) + 1013904223U;
    return (s >> 8) * (1.F / (1U << 24));
}

static float scnoise_val(Vec3 p, float cell_size, uint32_t seed) {
    const float inv_cs   = 1.F / cell_size;
    const float cs2      = cell_size * cell_size;
    const float inv_l2   = 9.F * inv_cs * inv_cs;
    const float inv_norm = 1.F / sqrtf(SC_KERNEL_VAR);

    const int gx = static_cast<int>(floorf(p.x * inv_cs));
    const int gy = static_cast<int>(floorf(p.y * inv_cs));
    const int gz = static_cast<int>(floorf(p.z * inv_cs));

    float vsum = 0.F;

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                const int cx = gx + dx;
                const int cy = gy + dy;
                const int cz = gz + dz;
                uint32_t cs = cell_hash(cx, cy, cz, seed);

                for (int k = 0; k < SC_IMPULSES; ++k) {
                    const float xi = cell_size * (static_cast<float>(cx) + lcg_next(cs));
                    const float yi = cell_size * (static_cast<float>(cy) + lcg_next(cs));
                    const float zi = cell_size * (static_cast<float>(cz) + lcg_next(cs));
                    const float wi = (lcg_next(cs) < 0.5F) ? 1.F : -1.F;

                    const float rx = p.x - xi;
                    const float ry = p.y - yi;
                    const float rz = p.z - zi;
                    const float d2 = (rx * rx) + (ry * ry) + (rz * rz);
                    if (d2 >= cs2) {
                        continue;
}
                    vsum += wi * sc_kernel(d2, inv_l2);
                }
            }
}
}
    return vsum * inv_norm;
}

static void
scnoise_val_grad(Vec3 p, float cell_size, uint32_t seed, float &out_val, Vec3 &out_grad) {
    const float inv_cs   = 1.F / cell_size;
    const float cs2      = cell_size * cell_size;
    const float inv_l2   = 9.F * inv_cs * inv_cs;
    const float inv_norm = 1.F / sqrtf(SC_KERNEL_VAR);

    const int gx = static_cast<int>(floorf(p.x * inv_cs));
    const int gy = static_cast<int>(floorf(p.y * inv_cs));
    const int gz = static_cast<int>(floorf(p.z * inv_cs));

    float vsum = 0.F;
    Vec3 gsum  = {0.F, 0.F, 0.F};

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                const int cx = gx + dx;
                const int cy = gy + dy;
                const int cz = gz + dz;
                uint32_t cs = cell_hash(cx, cy, cz, seed);

                for (int k = 0; k < SC_IMPULSES; ++k) {
                    const float xi = cell_size * (static_cast<float>(cx) + lcg_next(cs));
                    const float yi = cell_size * (static_cast<float>(cy) + lcg_next(cs));
                    const float zi = cell_size * (static_cast<float>(cz) + lcg_next(cs));
                    const float wi = (lcg_next(cs) < 0.5F) ? 1.F : -1.F;

                    const float rx = p.x - xi;
                    const float ry = p.y - yi;
                    const float rz = p.z - zi;
                    const float d2 = (rx * rx) + (ry * ry) + (rz * rz);
                    if (d2 >= cs2) {
                        continue;
}

                    float kv  = sc_kernel(d2, inv_l2);
                    float dkc = -inv_l2 * kv;
                    vsum += wi * kv;
                    gsum.x += wi * dkc * rx;
                    gsum.y += wi * dkc * ry;
                    gsum.z += wi * dkc * rz;
                }
            }
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
    float tx0 = (b.lower_x - ox) * idx;
    float tx1 = (b.upper_x - ox) * idx;
    float ty0 = (b.lower_y - oy) * idy;
    float ty1 = (b.upper_y - oy) * idy;
    float tz0 = (b.lower_z - oz) * idz;
    float tz1 = (b.upper_z - oz) * idz;
    if (tx0 > tx1) {
        std::swap(tx0, tx1);
}
    if (ty0 > ty1) {
        std::swap(ty0, ty1);
}
    if (tz0 > tz1) {
        std::swap(tz0, tz1);
}
    t0 = std::max({tx0, ty0, tz0, tnear});
    t1 = std::min({tx1, ty1, tz1, tfar});
    return t0 <= t1;
}

static RTCBounds segment_aabb(Vec3 p0, Vec3 p1, Vec3 p2, Vec3 p3, float r, float amplitude) {
    const float rr = r * (1.F + (3.F * amplitude));
    RTCBounds b;
    b.lower_x = std::min({p0.x, p1.x, p2.x, p3.x}) - rr;
    b.lower_y = std::min({p0.y, p1.y, p2.y, p3.y}) - rr;
    b.lower_z = std::min({p0.z, p1.z, p2.z, p3.z}) - rr;
    b.upper_x = std::max({p0.x, p1.x, p2.x, p3.x}) + rr;
    b.upper_y = std::max({p0.y, p1.y, p2.y, p3.y}) + rr;
    b.upper_z = std::max({p0.z, p1.z, p2.z, p3.z}) + rr;
    return b;
}

thread_local CondState tl_cond = {.u_grad={0.F, 0.F, 0.F}, .u_val=0.F, .valid=false};

struct MarchResult {
    float t_hit;
    float ct;
    Vec3 normal;
    float u_val;
    Vec3 u_grad;
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
                                 Vec3 p0,
                                 Vec3 p1,
                                 Vec3 p2,
                                 Vec3 p3,
                                 float r,
                                 const RTCBounds &aabb,
                                 float tnn,
                                 float amplitude,
                                 float cell_size,
                                 uint32_t seed) {
    MarchResult result{};

    float idx = (fabsf(dx) > 1e-9F) ? 1.F / dx : 1e30F * (dx >= 0.F ? 1.F : -1.F);
    float idy = (fabsf(dy) > 1e-9F) ? 1.F / dy : 1e30F * (dy >= 0.F ? 1.F : -1.F);
    float idz = (fabsf(dz) > 1e-9F) ? 1.F / dz : 1e30F * (dz >= 0.F ? 1.F : -1.F);

    float t0;
    float t1;
    if (!ray_aabb(ox, oy, oz, idx, idy, idz, tnn, tff, aabb, t0, t1)) {
        return result;
}

    const Vec3 o = {ox, oy, oz};
    const Vec3 d = {dx, dy, dz};

    constexpr int N = 16;
    Vec3 caps[N + 1];
    for (int i = 0; i <= N; ++i) {
        caps[i] = bspline_pos(p0, p1, p2, p3, i / static_cast<float>(N));
}

    const float l_noise = (amplitude > 0.F && cell_size > 1e-8F)
                        ? SC_IMPULSES * 27.F * 0.60653F * (3.F / cell_size) / sqrtf(SC_KERNEL_VAR)
                        : 0.F;
    const float l_total = 1.F + (amplitude * l_noise);

    const float hit_eps  = r * 1e-3F;
    const float min_step = r * 1e-5F;
    constexpr int MAX_STEPS = 256;

    float t = t0;

    for (int step = 0; step < MAX_STEPS; ++step) {
        if (t > t1) { break;
}

        const Vec3 q = o + d * t;

        float min_dist_sq = 1e30F;
        Vec3  closest_v   = {};

        for (int i = 0; i < N; ++i) {
            const Vec3  a    = caps[i];
            const Vec3  ab   = caps[i + 1] - a;
            const float len2 = dot(ab, ab);
            const float s    = (len2 > 1e-12F)
                             ? std::clamp(dot(q - a, ab) / len2, 0.F, 1.F)
                             : 0.F;
            const Vec3  v    = q - (a + ab * s);
            const float d2   = dot(v, v);
            if (d2 < min_dist_sq) { min_dist_sq = d2; closest_v = v; }
        }

        const float min_dist = sqrtf(min_dist_sq);

        float noise_val  = 0.F;
        Vec3  noise_grad = {};
        if (amplitude > 0.F) {
            scnoise_val_grad(q, cell_size, seed, noise_val, noise_grad);
}

        const float phi = min_dist - r - (amplitude * noise_val);

        if (phi <= hit_eps) {
            const Vec3 n_cap  = (min_dist > 1e-8F)
                              ? closest_v * (1.F / min_dist)
                              : Vec3{0.F, 1.F, 0.F};
            const Vec3 normal = normalize(n_cap - noise_grad * amplitude);
            return {.t_hit=t, .ct=0.F, .normal=normal, .u_val=0.F, .u_grad={0.F, 0.F, 0.F}, .gz=0.F, .noise_val=0.F, .hit=true};
        }

        t += std::max(phi / l_total, min_step);
    }

    return result;
}

static void bounds_cb(const RTCBoundsFunctionArguments *args) {
    const auto *hair = static_cast<const HairData *>(args->geometryUserPtr);
    int base             = hair->indices[args->primID];
    float r              = std::max({hair->vertices[base].w,
                                     hair->vertices[base + 1].w,
                                     hair->vertices[base + 2].w,
                                     hair->vertices[base + 3].w});
    Vec3 p0              = from_f4(hair->vertices[base + 0]);
    Vec3 p1              = from_f4(hair->vertices[base + 1]);
    Vec3 p2              = from_f4(hair->vertices[base + 2]);
    Vec3 p3              = from_f4(hair->vertices[base + 3]);
    *args->bounds_o      = segment_aabb(p0, p1, p2, p3, r, hair->amplitude);
}

static void intersect_cb(const RTCIntersectFunctionNArguments *args) {
    const auto *hair    = static_cast<const HairData *>(args->geometryUserPtr);
    const unsigned int n    = args->N;
    const unsigned int prim = args->primID;
    const unsigned int geom = args->geomID;

    RTCRayN *ray_n = RTCRayHitN_RayN(args->rayhit, n);
    RTCHitN *hit_n = RTCRayHitN_HitN(args->rayhit, n);

    int base = hair->indices[prim];
    float r  = std::max({hair->vertices[base].w,
                         hair->vertices[base + 1].w,
                         hair->vertices[base + 2].w,
                         hair->vertices[base + 3].w});
    Vec3 p0  = from_f4(hair->vertices[base + 0]);
    Vec3 p1  = from_f4(hair->vertices[base + 1]);
    Vec3 p2  = from_f4(hair->vertices[base + 2]);
    Vec3 p3  = from_f4(hair->vertices[base + 3]);

    RTCBounds aabb        = segment_aabb(p0, p1, p2, p3, r, hair->amplitude);
    const float amplitude = hair->amplitude * r;

    for (unsigned int i = 0; i < n; ++i) {
        if (!args->valid[i]) {
            continue;
}

        float ox  = RTCRayN_org_x(ray_n, n, i);
        float oy  = RTCRayN_org_y(ray_n, n, i);
        float oz  = RTCRayN_org_z(ray_n, n, i);
        float ddx = RTCRayN_dir_x(ray_n, n, i);
        float ddy = RTCRayN_dir_y(ray_n, n, i);
        float ddz = RTCRayN_dir_z(ray_n, n, i);
        float tnn = RTCRayN_tnear(ray_n, n, i);
        float tff = RTCRayN_tfar(ray_n, n, i);

        MarchResult res = march_segment(
            ox, oy, oz, ddx, ddy, ddz, tff, p0, p1, p2, p3, r, aabb, tnn,
            amplitude, hair->cell_size, hair->seed);
        if (!res.hit) {
            continue;
}

        RTCRayN_tfar(ray_n, n, i)   = res.t_hit;
        RTCHitN_Ng_x(hit_n, n, i)   = res.normal.x;
        RTCHitN_Ng_y(hit_n, n, i)   = res.normal.y;
        RTCHitN_Ng_z(hit_n, n, i)   = res.normal.z;
        RTCHitN_u(hit_n, n, i)      = res.ct;
        RTCHitN_v(hit_n, n, i)      = res.u_val;
        RTCHitN_primID(hit_n, n, i) = prim;
        RTCHitN_geomID(hit_n, n, i) = geom;
        tl_cond = {.u_grad=res.u_grad, .u_val=res.u_val, .valid=true};
    }
}

static void occluded_cb(const RTCOccludedFunctionNArguments *args) {
    const auto *hair    = static_cast<const HairData *>(args->geometryUserPtr);
    const unsigned int n    = args->N;
    const unsigned int prim = args->primID;

    RTCRayN *ray_n = args->ray;

    int base = hair->indices[prim];
    float r  = std::max({hair->vertices[base].w,
                         hair->vertices[base + 1].w,
                         hair->vertices[base + 2].w,
                         hair->vertices[base + 3].w});
    Vec3 p0  = from_f4(hair->vertices[base + 0]);
    Vec3 p1  = from_f4(hair->vertices[base + 1]);
    Vec3 p2  = from_f4(hair->vertices[base + 2]);
    Vec3 p3  = from_f4(hair->vertices[base + 3]);

    RTCBounds aabb        = segment_aabb(p0, p1, p2, p3, r, hair->amplitude);
    const float amplitude = hair->amplitude * r;

    for (unsigned int i = 0; i < n; ++i) {
        if (!args->valid[i]) {
            continue;
}

        float ox  = RTCRayN_org_x(ray_n, n, i);
        float oy  = RTCRayN_org_y(ray_n, n, i);
        float oz  = RTCRayN_org_z(ray_n, n, i);
        float ddx = RTCRayN_dir_x(ray_n, n, i);
        float ddy = RTCRayN_dir_y(ray_n, n, i);
        float ddz = RTCRayN_dir_z(ray_n, n, i);
        float tnn = RTCRayN_tnear(ray_n, n, i);
        float tff = RTCRayN_tfar(ray_n, n, i);

        MarchResult res = march_segment(
            ox, oy, oz, ddx, ddy, ddz, tff, p0, p1, p2, p3, r, aabb, tnn,
            amplitude, hair->cell_size, hair->seed);
        if (res.hit) {
            // bit-cast to avoid UB under -ffast-math
            static const uint32_t neg_inf_bits = 0xFF800000U;
            float neg_inf;
            __builtin_memcpy(&neg_inf, &neg_inf_bits, 4);
            RTCRayN_tfar(ray_n, n, i) = neg_inf;
        }
    }
}

GpisHitInfo extract_hit_info(const RTCRayHit &rh, const HairData &hair) {
    const float r  = hair.vertices[0].w;
    const float l2 = r * r / 9.F;

    return {.p={rh.ray.org_x + (rh.ray.dir_x * rh.ray.tfar),
             rh.ray.org_y + (rh.ray.dir_y * rh.ray.tfar),
             rh.ray.org_z + (rh.ray.dir_z * rh.ray.tfar)},
            .grad_mu=normalize({rh.hit.Ng_x, rh.hit.Ng_y, rh.hit.Ng_z}),
            .gz=0.F,
            .l2=l2,
            .kdd0=kappa_dd0(l2)};
}

unsigned add_user_hair(RTCDevice device, RTCScene scene, const HairData &hair) {
    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
    rtcSetGeometryUserPrimitiveCount(geom, static_cast<unsigned>(hair.indices.size()));
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
