#pragma once
#include <embree4/rtcore.h>
#include "hair_loader.h"

namespace m3hair {

// Attach a GPIS user geometry to an Embree scene for one HairData.
// Each B-spline segment becomes one user primitive.
// Returns the geomID assigned by Embree.
unsigned add_user_hair(RTCDevice device, RTCScene scene, const HairData& hair);

} // namespace m3hair
