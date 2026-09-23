# Standalone C ABI over the DspTap primitives (fft, yin, psola, pvoc, log_mel, decimate), used
# by the verification notebooks in ../../notebooks via ctypes. Builds on its own:
#
#     cmake -B build_capi -S tools/capi
#     cmake --build build_capi
#
# Copyright 2026 Timothy Place and the DspTap contributors. MIT License.

cmake_minimum_required(VERSION 3.19)
project(dsptap_capi CXX)

add_library(dsptap_capi SHARED dsptap_capi.cpp)

# The FFT entry points and pvoc pull the shared real FFT (tap::dsp). When this builds as part of the DspTap project the
# target already exists; standalone, pull the repo root in directly.
if (NOT TARGET tap::dsp)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/dsptap)
endif ()
target_link_libraries(dsptap_capi PRIVATE tap::dsp)
# dsptap_capi.h exports its symbols while this library builds and imports them everywhere else
# (the __declspec(dllexport) / dllimport split on Windows; visibility("default") elsewhere).
target_compile_definitions(dsptap_capi PRIVATE DSPTAP_BUILDING)
if (TARGET tap_dsp_warnings)
    target_link_libraries(dsptap_capi PRIVATE tap_dsp_warnings)
endif ()

set_target_properties(dsptap_capi PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    OUTPUT_NAME dsptap_capi
    CXX_VISIBILITY_PRESET hidden)

# The C-side smoke test (test_capi.c): compiles the header with a C front end and drives every
# FFT profile through the exponent contract the way a binding does. Built only when the root
# project builds its tests (never in the standalone build the ctypes bridge does), registered
# with CTest as dsptap_capi_smoke; it lives here rather than under tests/ so the executable lands
# next to the shared library on every platform (Windows needs the DLL beside the .exe).
if (TAP_DSP_BUILD_TESTS AND NOT CMAKE_CROSSCOMPILING)
    enable_language(C)
    add_executable(dsptap_capi_smoke test_capi.c)
    target_link_libraries(dsptap_capi_smoke PRIVATE dsptap_capi $<$<NOT:$<PLATFORM_ID:Windows>>:m>)
    set_target_properties(dsptap_capi_smoke PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
    add_test(NAME dsptap_capi_smoke COMMAND dsptap_capi_smoke)
endif ()
