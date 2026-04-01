#pragma once
#include "hair_loader.h"

#include <embree4/rtcore.h>

namespace m3hair {

class Scene {
  public:
    explicit Scene(RTCDevice device);
    ~Scene();

    void add_hair(const HairData &hair);
    void add_user_hair(const HairData &hair);
    void commit();

    RTCScene handle() const { return m_scene; }

  private:
    RTCDevice m_device;
    RTCScene m_scene;
};

} // namespace m3hair
