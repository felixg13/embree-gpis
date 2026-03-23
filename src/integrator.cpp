#include "integrator.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace m3hair {

// ---------------------------------------------------------------------------
// Cubic B-spline tangent at parameter u in [0,1] for segment with CPs p0..p3.
// Derivative of the de Boor basis: dB/du at u.
// ---------------------------------------------------------------------------
static vec3 bspline_tangent_at(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float u)
{
    // Derivative of cubic uniform B-spline basis
    // dB0/du = -u²/2 + u - 1/2
    // dB1/du = 3u²/2 - 2u
    // dB2/du = -3u²/2 + u + 1/2
    // dB3/du = u²/2
    float u2 = u*u;
    float d0 = -u2*0.5f + u       - 0.5f;
    float d1 =  u2*1.5f - 2.f*u;
    float d2 = -u2*1.5f + u       + 0.5f;
    float d3 =  u2*0.5f;
    return p0*d0 + p1*d1 + p2*d2 + p3*d3;
}

vec3 bspline_tangent(const HairData& hair, unsigned primID, float u)
{
    int base = hair.indices[primID];
    vec3 p0 = { hair.vertices[base+0].x, hair.vertices[base+0].y, hair.vertices[base+0].z };
    vec3 p1 = { hair.vertices[base+1].x, hair.vertices[base+1].y, hair.vertices[base+1].z };
    vec3 p2 = { hair.vertices[base+2].x, hair.vertices[base+2].y, hair.vertices[base+2].z };
    vec3 p3 = { hair.vertices[base+3].x, hair.vertices[base+3].y, hair.vertices[base+3].z };
    vec3 T = bspline_tangent_at(p0, p1, p2, p3, u);
    float len = length(T);
    return len > 1e-6f ? T / len : vec3{0.f, 1.f, 0.f};
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static RTCRayHit make_ray(const Ray& r) {
    RTCRayHit rh{};
    rh.ray.org_x = r.origin.x; rh.ray.org_y = r.origin.y; rh.ray.org_z = r.origin.z;
    rh.ray.dir_x = r.dir.x;   rh.ray.dir_y = r.dir.y;   rh.ray.dir_z = r.dir.z;
    rh.ray.tnear = r.tnear;   rh.ray.tfar  = r.tfar;
    rh.ray.mask  = -1;
    rh.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    return rh;
}

static float lum(const vec3& v) {
    return 0.2126f*v.x + 0.7152f*v.y + 0.0722f*v.z;
}

// ---------------------------------------------------------------------------
// Path tracer
// ---------------------------------------------------------------------------

vec3 trace_path(const Ray& ray,
                RTCScene scene,
                const DirectionalLight& light,
                const std::vector<const HairData*>& hairs,
                const std::vector<DeonParams>& params,
                int max_depth,
                RNG& rng)
{
    RTCIntersectArguments iargs; rtcInitIntersectArguments(&iargs);
    RTCOccludedArguments  oargs; rtcInitOccludedArguments(&oargs);

    vec3 L(0.f);
    vec3 throughput(1.f);
    Ray  cur = ray;

    for (int depth = 0; depth < max_depth; ++depth) {
        RTCRayHit rh = make_ray(cur);
        rtcIntersect1(scene, &rh, &iargs);

        if (rh.hit.geomID == RTC_INVALID_GEOMETRY_ID) break;   // miss → black bg

        // Select geometry data for this hit
        unsigned geom_idx = std::min((unsigned)rh.hit.geomID, (unsigned)(hairs.size()-1));
        const HairData& hair  = *hairs[geom_idx];
        const DeonParams& dp  = params[std::min(geom_idx, (unsigned)(params.size()-1))];

        vec3 T = bspline_tangent(hair, rh.hit.primID, rh.hit.u);

        vec3 hit_pos = {
            cur.origin.x + rh.ray.tfar * cur.dir.x,
            cur.origin.y + rh.ray.tfar * cur.dir.y,
            cur.origin.z + rh.ray.tfar * cur.dir.z
        };

        vec3 wo = -cur.dir;   // view direction (toward camera)

        // --- NEE: shadow ray toward light ---
        {
            RTCRay shadow{};
            shadow.org_x = hit_pos.x; shadow.org_y = hit_pos.y; shadow.org_z = hit_pos.z;
            shadow.dir_x = light.direction.x;
            shadow.dir_y = light.direction.y;
            shadow.dir_z = light.direction.z;
            shadow.tnear = 1e-3f;
            shadow.tfar  = 1e30f;
            shadow.mask  = -1;

            rtcOccluded1(scene, &shadow, &oargs);

            if (shadow.tfar > 0.f) {   // not occluded (Embree sets tfar=-inf on hit)
                vec3 wi_light = light.direction;
                vec3 f = eval_deon(wi_light, wo, T, dp);
                float cos_o = std::abs(dot(wi_light, T));
                L = L + throughput * f * light.radiance * cos_o;
            }
        }

        // --- Russian roulette at depth >= 3 ---
        if (depth >= 3) {
            float q = std::min(lum(throughput), 0.95f);
            if (rng.next_f() > q) break;
            throughput = throughput * (1.f / q);
        }

        // --- Sample BSDF for next ray ---
        vec3 wi_new;
        float pdf;
        vec3 f = sample_deon(wo, T, dp,
                             rng.next_f(), rng.next_f(), rng.next_f(), rng.next_f(),
                             wi_new, pdf);

        if (pdf < 1e-7f) break;

        float cos_i = std::abs(dot(wi_new, T));
        throughput = throughput * f * (cos_i / pdf);

        if (lum(throughput) < 1e-4f) break;

        cur = Ray{ hit_pos, wi_new, 1e-3f, 1e30f };
    }

    return L;
}

} // namespace m3hair
