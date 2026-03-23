cmake_minimum_required(VERSION 3.20)

project(gpis_hair
    VERSION 0.1.0
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-O3 -march=native -Wall -Wextra -ffast-math)
endif()

find_package(embree 4 REQUIRED)
find_package(TBB REQUIRED)

# Shared library sources (no main)
set(LIB_SOURCES
    src/hair_loader.cpp
    src/scene.cpp
    src/image.cpp
    src/bsdf_deon.cpp
    src/integrator.cpp
    src/user_geo.cpp
)

# Main renderer
add_executable(gpis_hair src/main.cpp ${LIB_SOURCES})
target_include_directories(gpis_hair PRIVATE src)
target_link_libraries(gpis_hair PRIVATE embree4 TBB::tbb)

# BSDF unit tests
add_executable(bsdf_test src/bsdf_test.cpp ${LIB_SOURCES})
target_include_directories(bsdf_test PRIVATE src)
target_link_libraries(bsdf_test PRIVATE embree4 TBB::tbb)
