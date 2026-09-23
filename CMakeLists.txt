# 3.22 is Ubuntu 22.04 LTS's stock cmake. Nothing here uses a 3.23+ feature,
# so the higher floor only excluded the current LTS.
cmake_minimum_required(VERSION 3.22)
project(DspTap VERSION 0.1.0 LANGUAGES C CXX)

# ==============================================================================
# DspTap — shared DSP primitives for the Tap family of audio libraries.
#
# Header-only C++20 (include/tap/dsp/): the real FFT in four profiles — the
# C++20 split-radix engine for double and float (fft/split_radix.h, the
# transliteration of Ooura's rdft that replaced the vendored C at Stage 2b and
# is bit-identical to it), the int32 fixed-point kernel for Q15 and Q31
# (fft/fixed_point.h) — with optional, mutually-exclusive float32 backends
# (CMSIS-DSP Helium on the Cortex-M55, Apple vDSP on macOS) that re-present
# the same numeric contract; plus the pitch, wake-word, inference and FIR
# primitives README.md lists. Consumed as a submodule by MuTap, AmbiTap, ….
#
# tap::dsp is a pure INTERFACE target. Since Stage 2c
# (docs/audit-fft-and-code-smells.md, Part 3; Decision D6) no vendored C is
# compiled into what ships: the reference copy of Ooura's fftsg.c lives under
# tests/reference/ooura/ and is compiled only by the parity gate in tests/.
# The one compiled object library, tap_dsp_fft, exists only when
# TAP_DSP_FFT_CMSIS is ON and carries the CMSIS-DSP objects (below).
# ==============================================================================

# ------------------------------------------------------------------------------
# The library proper (header-only)
# ------------------------------------------------------------------------------

add_library(tap_dsp INTERFACE)
add_library(tap::dsp ALIAS tap_dsp)
target_include_directories(tap_dsp INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(tap_dsp INTERFACE cxx_std_20)

# ------------------------------------------------------------------------------
# Helium (Cortex-M55 MVE) FFT backend. When ON, the float32 real FFT in
# include/tap/dsp/fft.h routes through the vendored CMSIS-DSP subset
# (third_party/cmsis-dsp, Apache-2.0) instead of the split-radix engine —
# measured 1.8-2.0x fewer instructions over the whole ratchet scenario on the
# M55 (bench/README.md). Double is untouched (always the split-radix engine),
# so the double-precision golden model is unaffected; on the float32 embedded
# path CMSIS agrees with the engine to single-precision rounding (not
# bit-identical). This is the one case in which tap::dsp links a compiled
# library: tap_dsp_fft (alias tap::dsp_fft), PIC so MODULE/SHARED consumers
# (host applications, Max externals) can link it, carrying the CMSIS objects
# and nothing else. Until Stage 2c the same target also carried the vendored
# Ooura C on every platform; it does not exist at all without this option now.
#
# Default ON for the bare-metal M55 embedded profile (cmake/arm-cortex-m55-mps3.cmake),
# OFF everywhere else — desktop, the Max/C-ABI host builds (incl. Apple Silicon
# arm64, which is not this Generic+arm profile), and Hexagon all keep the
# split-radix engine. Force it on the M55 with -DTAP_DSP_FFT_CMSIS=OFF.
# ------------------------------------------------------------------------------
if(CMAKE_SYSTEM_NAME STREQUAL "Generic" AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm|ARM")
    set(_tap_dsp_fft_cmsis_default ON) # bare-metal M55 (cmake/arm-cortex-m55-mps3.cmake)
else()
    set(_tap_dsp_fft_cmsis_default OFF)
endif()
option(TAP_DSP_FFT_CMSIS "Route the float32 FFT through the vendored CMSIS-DSP Helium backend (Arm/Helium only)" ${_tap_dsp_fft_cmsis_default})
if(TAP_DSP_FFT_CMSIS)
    if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "arm|ARM")
        message(FATAL_ERROR
            "TAP_DSP_FFT_CMSIS is an Arm/Helium backend; configure with the "
            "cortex-m55 toolchain (cmake/arm-cortex-m55-mps3.cmake).")
    endif()
    file(GLOB _tap_dsp_cmsis_src CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/Source/*/*.c)
    add_library(tap_dsp_fft STATIC ${_tap_dsp_cmsis_src})
    add_library(tap::dsp_fft ALIAS tap_dsp_fft)
    set_target_properties(tap_dsp_fft PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        OUTPUT_NAME tap_dsp_fft)
    target_link_libraries(tap_dsp INTERFACE tap_dsp_fft)
    # PRIVATE for the vendored .c compiles; the same dirs are exposed on the
    # INTERFACE below so consumers' fft.h can resolve "arm_math.h".
    set(_tap_dsp_cmsis_inc
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/Include
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/PrivateInclude
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/cmsis-core/Include)
    target_include_directories(tap_dsp_fft PRIVATE ${_tap_dsp_cmsis_inc})
    target_include_directories(tap_dsp INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/Include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/PrivateInclude>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/cmsis-core/Include>)
    target_compile_definitions(tap_dsp INTERFACE TAP_DSP_FFT_CMSIS)
    message(STATUS "tap::dsp: float32 FFT backend = CMSIS-DSP Helium")
endif()

# ------------------------------------------------------------------------------
# Apple vDSP (Accelerate) FFT backend. When ON, the float32 real FFT in
# include/tap/dsp/fft.h routes through Apple's vDSP instead of the split-radix
# engine — measured ~3x faster per transform than the clang-autovectorized C
# on Apple Silicon (MuTap's transform-only figure; the same-binary comparison
# is Stage 4). Double is untouched (always the split-radix engine), so the
# double-precision golden model is unaffected; on the float32 path vDSP agrees
# with the engine to single-precision rounding (not bit-identical). Default ON
# for Apple, OFF everywhere else; it is mutually exclusive with the CMSIS
# backend. No compiled library: the framework is linked on the INTERFACE.
# Force the split-radix engine with -DTAP_DSP_FFT_ACCELERATE=OFF.
# ------------------------------------------------------------------------------
if(APPLE AND NOT TAP_DSP_FFT_CMSIS)
    set(_tap_dsp_fft_accel_default ON)
else()
    set(_tap_dsp_fft_accel_default OFF)
endif()
option(TAP_DSP_FFT_ACCELERATE "Route the float32 FFT through Apple's vDSP (Accelerate) backend (Apple only)" ${_tap_dsp_fft_accel_default})
if(TAP_DSP_FFT_ACCELERATE)
    if(NOT APPLE)
        message(FATAL_ERROR "TAP_DSP_FFT_ACCELERATE requires Apple's Accelerate framework (macOS).")
    endif()
    if(TAP_DSP_FFT_CMSIS)
        message(FATAL_ERROR "TAP_DSP_FFT_ACCELERATE and TAP_DSP_FFT_CMSIS are mutually exclusive.")
    endif()
    target_compile_definitions(tap_dsp INTERFACE TAP_DSP_FFT_ACCELERATE)
    target_link_libraries(tap_dsp INTERFACE "-framework Accelerate")
    message(STATUS "tap::dsp: float32 FFT backend = Apple vDSP (Accelerate)")
endif()

# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------

if(PROJECT_IS_TOP_LEVEL)
    option(TAP_DSP_BUILD_TESTS "Build DspTap tests" ON)
else()
    option(TAP_DSP_BUILD_TESTS "Build DspTap tests" OFF)
endif()

# Warning flags for this project's own tests; never exported to consumers of
# the INTERFACE library, and never applied to vendored C (the CMSIS objects
# under TAP_DSP_FFT_CMSIS, the reference C the parity gate compiles).
option(TAP_DSP_WERROR "Treat warnings as errors in DspTap's own targets" OFF)
add_library(tap_dsp_warnings INTERFACE)
target_compile_options(tap_dsp_warnings INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic -Wconversion -Wshadow>
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive->)
if(TAP_DSP_WERROR)
    target_compile_options(tap_dsp_warnings INTERFACE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Werror>
        $<$<CXX_COMPILER_ID:MSVC>:/WX>)
endif()

if(TAP_DSP_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Benchmarks: the host wall-clock microbenchmark (informational) and the
# deterministic instruction-count ratchet workloads (bench/icount), one binary
# per scenario, counted under QEMU by scripts/icount.py from
# .github/workflows/bench.yml. Builds on the host and on the cross toolchains;
# policy and workflow in bench/README.md.
option(TAP_DSP_BUILD_BENCH "Build DspTap benchmarks and instruction-count ratchet workloads" OFF)
if(TAP_DSP_BUILD_BENCH)
    add_subdirectory(bench)
endif()

# C ABI shared library for FFI consumers (the notebooks drive the shipping C++
# through it via ctypes — see notebooks/dsptap_py.py). tools/capi also builds
# standalone (`cmake -B build_capi -S tools/capi`), which is what the ctypes
# bridge does on first import; this option exists so the root build can compile
# it too, matching the BUILD_CAPI option every sibling library carries. CI turns
# it on so the verification layer cannot rot unnoticed.
option(TAP_DSP_BUILD_CAPI "Build the C ABI shared library" OFF)
if(TAP_DSP_BUILD_CAPI)
    add_subdirectory(tools/capi)
endif()
