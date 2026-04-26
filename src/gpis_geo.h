#pragma once
#include "gpis_nee.h"
#include "hair_loader.h"
#include "math.h"

#include <embree4/rtcore.h>

namespace m3hair {

struct CondState {
    Vec3 u_grad;
    float u_val;
    bool valid;
};
extern thread_local CondState tl_cond;

unsigned add_user_hair(RTCDevice device, RTCScene scene, const HairData &hair);
GpisHitInfo extract_hit_info(const RTCRayHit &rh, const HairData &hair);

} // namespace m3hair
