cmake_minimum_required(VERSION 3.20)

project(gpis_hair
    VERSION 0.1.0
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(CMAKE_BUILD_TYPE STREQUAL "Profile")
        add_compile_options(-O2 -g -fno-omit-frame-pointer -march=native -Wall -Wextra)
        add_compile_definitions(ENABLE_PERF_COUNTERS)
    else()
        add_compile_options(-O3 -march=native -Wall -Wextra -ffast-math -fno-omit-frame-pointer)
    endif()
endif()

find_package(embree 4 REQUIRED)
find_package(TBB REQUIRED)

include(FetchContent)
FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG v1.16.0
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(spdlog)

# spdlog's bundled fmt uses isfinite/isnan which triggers -Wnan-infinity-disabled
# under -ffast-math.  Mark its headers as system includes to suppress the noise.
get_target_property(_spdlog_incs spdlog_header_only INTERFACE_INCLUDE_DIRECTORIES)
set_target_properties(spdlog_header_only PROPERTIES
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_spdlog_incs}")

# Shared library sources (no main)
set(LIB_SOURCES
    src/hair_loader.cpp
    src/scene.cpp
    src/image.cpp
    src/bsdf_deon.cpp
    src/integrator.cpp
    src/user_geo.cpp
    src/gpis_nee.cpp
)

# Main renderer
add_executable(gpis_hair src/main.cpp ${LIB_SOURCES})
target_include_directories(gpis_hair PRIVATE src)
target_link_libraries(gpis_hair PRIVATE embree4 TBB::tbb spdlog::spdlog_header_only)

# BSDF unit tests
add_executable(bsdf_test src/bsdf_test.cpp ${LIB_SOURCES})
target_include_directories(bsdf_test PRIVATE src)
target_link_libraries(bsdf_test PRIVATE embree4 TBB::tbb spdlog::spdlog_header_only)
