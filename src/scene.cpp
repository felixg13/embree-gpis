#include "scene.h"

#include "gpis_geo.h"

#include <stdexcept>

namespace m3hair {

Scene::Scene(RTCDevice device) : m_device(device) {
    m_scene = rtcNewScene(m_device);
    if (!m_scene)
        throw std::runtime_error("rtcNewScene failed");
}

Scene::~Scene() {
    rtcReleaseScene(m_scene);
}

void Scene::add_hair(const HairData &hair) {
    RTCGeometry geom = rtcNewGeometry(m_device, RTC_GEOMETRY_TYPE_ROUND_BSPLINE_CURVE);

    rtcSetSharedGeometryBuffer(geom,
                               RTC_BUFFER_TYPE_VERTEX,
                               0,
                               RTC_FORMAT_FLOAT4,
                               hair.vertices.data(),
                               0,
                               sizeof(float4),
                               hair.vertices.size());

    rtcSetSharedGeometryBuffer(geom,
                               RTC_BUFFER_TYPE_INDEX,
                               0,
                               RTC_FORMAT_UINT,
                               hair.indices.data(),
                               0,
                               sizeof(int),
                               hair.indices.size());

    rtcCommitGeometry(geom);
    rtcAttachGeometry(m_scene, geom);
    rtcReleaseGeometry(geom);
}

void Scene::add_user_hair(const HairData &hair) {
    m3hair::add_user_hair(m_device, m_scene, hair);
}

void Scene::commit() {
    rtcCommitScene(m_scene);
}

} // namespace m3hair
