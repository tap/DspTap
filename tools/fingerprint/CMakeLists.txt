# Same-host A/B fingerprint tool for the DspTap primitives whose output must not
# move when their internals do (pvoc, log_mel). Prints one FNV-1a-64 line per
# case; build it against two trees on ONE host with ONE compiler and flag set,
# and diff the two outputs. The hashes are host + compiler + flags specific
# (pvoc goes through atan2/sin/cos in libm; contraction changes them too), so
# they are never committed as golden values and never asserted in a test.
# Builds on its own:
#
#     cmake -B build_fp -S tools/fingerprint -DCMAKE_BUILD_TYPE=Release
#     cmake --build build_fp && build_fp/dsptap_fingerprint > after.txt
#
# Copyright 2026 Timothy Place and the DspTap contributors. MIT License.

cmake_minimum_required(VERSION 3.22)
project(dsptap_fingerprint CXX)

add_executable(dsptap_fingerprint fingerprint.cpp)

# pvoc and log_mel pull the shared real FFT (tap::dsp). When this builds as part
# of the DspTap project the target already exists; standalone, pull the repo
# root in directly (its tests stay off: it is not the top-level project).
if (NOT TARGET tap::dsp)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/dsptap)
endif ()
target_link_libraries(dsptap_fingerprint PRIVATE tap::dsp)
if (TARGET tap_dsp_warnings)
    target_link_libraries(dsptap_fingerprint PRIVATE tap_dsp_warnings)
endif ()

set_target_properties(dsptap_fingerprint PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF)
