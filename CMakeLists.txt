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
# compiled into what ships; the reference copy of Ooura's fftsg.c that the
# parity gate in tests/ compiled was deleted at D6, and third_party/ooura/
# readme.txt stays as the license record for the derived engine.
# The one compiled (static) library, tap_dsp_fft, exists only when
# TAP_DSP_FFT_CMSIS is ON and carries the CMSIS-DSP objects (below).
# ==============================================================================

# ------------------------------------------------------------------------------
# The library proper (header-only)
# ------------------------------------------------------------------------------

# Consumed by add_subdirectory only: there is no install() / export() rule, by
# decision (docs/audit-fft-and-code-smells.md, Part 6, P11), so the include
# directory has no $<INSTALL_INTERFACE> half (Stage 6 removed the dead one).
add_library(tap_dsp INTERFACE)
add_library(tap::dsp ALIAS tap_dsp)
target_include_directories(tap_dsp INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
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
# Default: ON exactly when the target compiler provides floating-point Helium
# (MVE-F) — the Cortex-M55 / M85 class — and OFF everywhere else: desktop, the
# Max/C-ABI host builds (incl. Apple Silicon arm64), Hexagon, and every Cortex-M
# without MVE-F (M0+/M4/M7/M33). Detected by a compile check on the compiler's
# own predefined macro, __ARM_FEATURE_MVE with the floating-point bit (2) set,
# under the configured flags (CMAKE_CXX_FLAGS, i.e. the toolchain file's
# CMAKE_CXX_FLAGS_INIT) — so it does not matter how a consumer's toolchain
# spells the CPU (-mcpu=cortex-m55, -mcpu=cortex-m85,
# -march=armv8.1-m.main+mve.fp, a CMAKE_SYSTEM_PROCESSOR of "arm" or anything
# else). Measured with arm-none-eabi-gcc 13.2: -mcpu=cortex-m55/m85 with
# -mfloat-abi=hard define it as 3; cortex-m55+nomve.fp as 1 (integer MVE only,
# default OFF); cortex-m55 with -mfloat-abi=soft, cortex-m33, cortex-m4(f) and
# cortex-m7 not at all. Until this check (audit build.md F3) the default was ON
# for every CMAKE_SYSTEM_NAME Generic + arm target, so a consumer's own M7 or
# M33 toolchain silently got the CMSIS engine: float size range 32 … 4096
# instead of 4 … 2^30, and the fft_cmsis ABI tag.
#
# Flags the check cannot see (per-target or per-directory options added after
# the toolchain) leave the default OFF: the split-radix engine is correct on
# every target, so a missed detection costs speed, never the contract. An
# explicit -DTAP_DSP_FFT_CMSIS=ON / OFF always wins. ON without MVE-F detected
# is allowed on a 32-bit Arm target, with a warning (CMSIS-DSP then compiles
# its non-Helium C, same 32 … 4096 range); ON for anything that is not 32-bit
# Arm (x86, arm64/Apple Silicon) is a configure error.
# ------------------------------------------------------------------------------
include(CheckCXXSourceCompiles)
# A compile-only check: a static-library try_compile never links, so it works
# on a bare-metal toolchain whether or not the toolchain file set
# CMAKE_TRY_COMPILE_TARGET_TYPE (no startup files or linker script needed).
# Function scope keeps both settings away from the caller's try_compiles.
function(_tap_dsp_target_defines condition result)
    set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
    set(CMAKE_REQUIRED_QUIET ON)
    check_cxx_source_compiles("#if !(${condition})\n#error \"${result}: no\"\n#endif\nint main() { return 0; }\n" ${result})
endfunction()
_tap_dsp_target_defines("defined(__ARM_FEATURE_MVE) && (__ARM_FEATURE_MVE & 2)" TAP_DSP_HAVE_MVE_FLOAT)
if(TAP_DSP_HAVE_MVE_FLOAT)
    set(_tap_dsp_fft_cmsis_default ON) # MVE-F: Cortex-M55 (cmake/arm-cortex-m55-mps3.cmake), M85
else()
    set(_tap_dsp_fft_cmsis_default OFF)
endif()
option(TAP_DSP_FFT_CMSIS "Route the float32 FFT through the vendored CMSIS-DSP Helium backend (default: ON when the compiler targets MVE-F)" ${_tap_dsp_fft_cmsis_default})
if(TAP_DSP_FFT_CMSIS)
    if(NOT TAP_DSP_HAVE_MVE_FLOAT)
        _tap_dsp_target_defines("defined(__arm__)" TAP_DSP_TARGET_IS_ARM32)
        if(NOT TAP_DSP_TARGET_IS_ARM32)
            message(FATAL_ERROR
                "TAP_DSP_FFT_CMSIS is a 32-bit Arm (Helium) backend and this compiler "
                "does not target 32-bit Arm; configure with the cortex-m55 toolchain "
                "(cmake/arm-cortex-m55-mps3.cmake).")
        endif()
        message(WARNING
            "TAP_DSP_FFT_CMSIS=ON but the compiler does not target floating-point "
            "Helium (__ARM_FEATURE_MVE & 2): CMSIS-DSP compiles its non-Helium C and "
            "the float FFT's size range is 32 … 4096 (split-radix: 4 … 2^30).")
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
    if(TAP_DSP_HAVE_MVE_FLOAT)
        message(STATUS "tap::dsp: float32 FFT backend = CMSIS-DSP Helium")
    else()
        message(STATUS "tap::dsp: float32 FFT backend = CMSIS-DSP (no Helium on this target)")
    endif()
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
# under TAP_DSP_FFT_CMSIS).
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
