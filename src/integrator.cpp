#include "integrator.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace m3hair {

// dB/du for cubic uniform B-spline
static Vec3 bspline_tangent_at(Vec3 p0, Vec3 p1, Vec3 p2, Vec3 p3, float u) {
    float u2 = u * u;
    float d0 = (-u2 * 0.5F) + u - 0.5F;
    float d1 = (u2 * 1.5F) - (2.F * u);
    float d2 = (-u2 * 1.5F) + u + 0.5F;
    float d3 = u2 * 0.5F;
    return p0 * d0 + p1 * d1 + p2 * d2 + p3 * d3;
}

Vec3 bspline_tangent(const HairData &hair, unsigned primID, float u) {
    int base  = hair.indices[primID];
    Vec3 p0   = {hair.vertices[base + 0].x, hair.vertices[base + 0].y, hair.vertices[base + 0].z};
    Vec3 p1   = {hair.vertices[base + 1].x, hair.vertices[base + 1].y, hair.vertices[base + 1].z};
    Vec3 p2   = {hair.vertices[base + 2].x, hair.vertices[base + 2].y, hair.vertices[base + 2].z};
    Vec3 p3   = {hair.vertices[base + 3].x, hair.vertices[base + 3].y, hair.vertices[base + 3].z};
    Vec3 t    = bspline_tangent_at(p0, p1, p2, p3, u);
    float len = length(t);
    return len > 1e-6F ? t / len : Vec3{0.F, 1.F, 0.F};
}

static RTCRayHit make_ray(const Ray &r) {
    RTCRayHit rh{};
    rh.ray.org_x  = r.origin.x;
    rh.ray.org_y  = r.origin.y;
    rh.ray.org_z  = r.origin.z;
    rh.ray.dir_x  = r.dir.x;
    rh.ray.dir_y  = r.dir.y;
    rh.ray.dir_z  = r.dir.z;
    rh.ray.tnear  = r.tnear;
    rh.ray.tfar   = r.tfar;
    rh.ray.mask   = -1;
    rh.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    return rh;
}

static float lum(const Vec3 &v) {
    return (0.2126F * v.x) + (0.7152F * v.y) + (0.0722F * v.z);
}

Vec3 trace_path(const Ray &ray,
                RTCScene scene,
                const DirectionalLight &light,
                const std::vector<const HairData *> &hairs,
                const std::vector<DeonParams> &params,
                int max_depth,
                RNG &rng) {
    RTCIntersectArguments iargs;
    rtcInitIntersectArguments(&iargs);
    RTCOccludedArguments oargs;
    rtcInitOccludedArguments(&oargs);

    Vec3 l(0.F);
    Vec3 throughput(1.F);
    Ray cur = ray;

    for (int depth = 0; depth < max_depth; ++depth) {
        RTCRayHit rh = make_ray(cur);
        rtcIntersect1(scene, &rh, &iargs);

        if (rh.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
            break;
}

        unsigned geom_idx    = std::min((unsigned)rh.hit.geomID, static_cast<unsigned>(hairs.size() - 1));
        const HairData &hair = *hairs[geom_idx];
        const DeonParams &dp = params[std::min(geom_idx, static_cast<unsigned>(params.size() - 1))];

        Vec3 t = bspline_tangent(hair, rh.hit.primID, rh.hit.u);

        Vec3 hit_pos = {cur.origin.x + (rh.ray.tfar * cur.dir.x),
                        cur.origin.y + (rh.ray.tfar * cur.dir.y),
                        cur.origin.z + (rh.ray.tfar * cur.dir.z)};

        Vec3 wo = -cur.dir;

        {
            RTCRay shadow{};
            shadow.org_x = hit_pos.x;
            shadow.org_y = hit_pos.y;
            shadow.org_z = hit_pos.z;
            shadow.dir_x = light.direction.x;
            shadow.dir_y = light.direction.y;
            shadow.dir_z = light.direction.z;
            shadow.tnear = 1e-3F;
            shadow.tfar  = 1e30F;
            shadow.mask  = -1;

            rtcOccluded1(scene, &shadow, &oargs);

            if (shadow.tfar > 0.F) {
                Vec3 f      = eval_deon(light.direction, wo, t, dp);
                float cos_o = std::abs(dot(light.direction, t));
                l           = l + throughput * f * light.radiance * cos_o;
            }
        }

        if (depth >= 3) {
            float q = std::min(lum(throughput), 0.95F);
            if (rng.next_f() > q) {
                break;
}
            throughput = throughput * (1.F / q);
        }

        Vec3 wi_new;
        float pdf;
        Vec3 f = sample_deon(
            wo, t, dp, rng.next_f(), rng.next_f(), rng.next_f(), rng.next_f(), wi_new, pdf);

        if (pdf < 1e-7F) {
            break;
}

        float cos_i = std::abs(dot(wi_new, t));
        throughput  = throughput * f * (cos_i / pdf);

        if (lum(throughput) < 1e-4F) {
            break;
}

        cur = Ray{.origin=hit_pos, .dir=wi_new, .tnear=1e-3F, .tfar=1e30F};
    }

    return l;
}

} // namespace m3hair
