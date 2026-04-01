#pragma once
#include "gpis_nee.h"
#include "hair_loader.h"
#include "math.h"

#include <embree4/rtcore.h>

namespace m3hair {

// ---------------------------------------------------------------------------
// Pathwise conditioning state.
// The renderer must set tl_cond before each rtcIntersect / rtcOccluded call:
//   - Camera rays (depth=0):    tl_cond = {{}, 0.f, false}
//   - Secondary rays (depth>0): tl_cond = {u_grad, u_val, true}
//     where u_grad and u_val were written by the previous intersect_cb.
// After rtcIntersect the callback writes new u_grad into tl_cond.u_grad and
// u_val into RTCHit.v so the renderer can pass them to the next bounce.
// ---------------------------------------------------------------------------
struct CondState {
    vec3 u_grad; // gradient conditioning weights
    float u_val; // value conditioning weight
    bool valid;  // false on camera rays (depth=0)
};
extern thread_local CondState tl_cond;

// ---------------------------------------------------------------------------
// Attach a GPIS user geometry to an Embree scene for one HairData.
// Each B-spline segment becomes one user primitive.
// Returns the geomID assigned by Embree.
// ---------------------------------------------------------------------------
unsigned add_user_hair(RTCDevice device, RTCScene scene, const HairData &hair);

// ---------------------------------------------------------------------------
// Reconstruct GpisHitInfo from an RTCRayHit record + the current HairData.
// Call after a successful rtcIntersect1 / rtcIntersectN.
// Note: info.gz is zeroed — fill it from your saved MarchResult or recompute.
// ---------------------------------------------------------------------------
GpisHitInfo extract_hit_info(const RTCRayHit &rh, const HairData &hair);

} // namespace m3hair
