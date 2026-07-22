# Standalone C ABI over the DspTap primitives (yin, psola, pvoc), used by the verification
# notebooks in ../../notebooks via ctypes. Builds on its own:
#
#     cmake -B build_capi -S tools/capi
#     cmake --build build_capi
#
# Copyright 2026 Timothy Place and the DspTap contributors. MIT License.

cmake_minimum_required(VERSION 3.19)
project(dsptap_capi CXX)

add_library(dsptap_capi SHARED dsptap_capi.cpp)

# pvoc pulls the shared real FFT (tap::dsp). When this builds as part of the DspTap project the
# target already exists; standalone, pull the repo root in directly.
if (NOT TARGET tap::dsp)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/dsptap)
endif ()
target_link_libraries(dsptap_capi PRIVATE tap::dsp)

set_target_properties(dsptap_capi PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    OUTPUT_NAME dsptap_capi
    CXX_VISIBILITY_PRESET hidden)
