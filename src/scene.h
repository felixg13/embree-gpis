#pragma once
#include <embree4/rtcore.h>
#include "hair_loader.h"

namespace m3hair {

class Scene {
public:
    explicit Scene(RTCDevice device);
    ~Scene();

    // Build (or rebuild) the Embree scene from hair data.
    // Multiple calls to add_hair() followed by commit().
    void add_hair(const HairData& hair);
    void commit();

    RTCScene handle() const { return m_scene; }

private:
    RTCDevice m_device;
    RTCScene  m_scene;
};

} // namespace m3hair
