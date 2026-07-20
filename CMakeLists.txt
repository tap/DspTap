cmake_minimum_required(VERSION 3.24)
project(DspTap VERSION 0.1.0 LANGUAGES C CXX)

# ==============================================================================
# DspTap — shared real-FFT primitive for the Tap family of DSP libraries.
#
# A header-only C++ wrapper (include/dsptap/fft.h) over the vendored Ooura
# split-radix real FFT, with optional, mutually-exclusive float32 backends
# (CMSIS-DSP Helium on the Cortex-M55, Apple vDSP on macOS) that re-present the
# same Ooura numeric contract. Consumed as a submodule by MuTap, AmbiTap, ….
# ==============================================================================

# ------------------------------------------------------------------------------
# Ooura FFT (vendored) — provides the extern "C" rdft/rdft_f used by
# include/dsptap/fft.h. PIC so MODULE/SHARED consumers (host applications,
# Max externals) can link it.
# ------------------------------------------------------------------------------

add_library(dsptap_fft STATIC third_party/ooura/fftsg.c third_party/ooura/fftsg_float.c)
add_library(DspTap::fft ALIAS dsptap_fft)
set_target_properties(dsptap_fft PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    OUTPUT_NAME dsptap_fft)

# ------------------------------------------------------------------------------
# The library proper (header-only)
# ------------------------------------------------------------------------------

add_library(dsptap INTERFACE)
add_library(DspTap::DspTap ALIAS dsptap)
target_include_directories(dsptap INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(dsptap INTERFACE cxx_std_20)
target_link_libraries(dsptap INTERFACE dsptap_fft)

# ------------------------------------------------------------------------------
# Helium (Cortex-M55 MVE) FFT backend. When ON, the float32 real FFT in
# include/dsptap/fft.h routes through the vendored CMSIS-DSP subset
# (third_party/cmsis-dsp, Apache-2.0) instead of Ooura — measured ~3x fewer
# instructions per transform on the M55. Double is untouched (always Ooura), so
# the double-precision golden model is unaffected; on the float32 embedded path
# CMSIS agrees with Ooura to single-precision rounding (not bit-identical).
#
# Default ON for the bare-metal M55 embedded profile (cmake/arm-cortex-m55-mps3.cmake),
# OFF everywhere else — desktop, the Max/C-ABI host builds (incl. Apple Silicon
# arm64, which is not this Generic+arm profile), and Hexagon all keep Ooura.
# Force the Ooura path on the M55 with -DDSPTAP_FFT_CMSIS=OFF.
# ------------------------------------------------------------------------------
if(CMAKE_SYSTEM_NAME STREQUAL "Generic" AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm|ARM")
    set(_dsptap_fft_cmsis_default ON) # bare-metal M55 (cmake/arm-cortex-m55-mps3.cmake)
else()
    set(_dsptap_fft_cmsis_default OFF)
endif()
option(DSPTAP_FFT_CMSIS "Route the float32 FFT through the vendored CMSIS-DSP Helium backend (Arm/Helium only)" ${_dsptap_fft_cmsis_default})
if(DSPTAP_FFT_CMSIS)
    if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "arm|ARM")
        message(FATAL_ERROR
            "DSPTAP_FFT_CMSIS is an Arm/Helium backend; configure with the "
            "cortex-m55 toolchain (cmake/arm-cortex-m55-mps3.cmake).")
    endif()
    file(GLOB _dsptap_cmsis_src CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/Source/*/*.c)
    target_sources(dsptap_fft PRIVATE ${_dsptap_cmsis_src})
    # PRIVATE for the vendored .c compiles; the same dirs are exposed on the
    # INTERFACE below so consumers' fft.h can resolve "arm_math.h".
    set(_dsptap_cmsis_inc
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/Include
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/PrivateInclude
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/cmsis-core/Include)
    target_include_directories(dsptap_fft PRIVATE ${_dsptap_cmsis_inc})
    target_include_directories(dsptap INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/Include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/PrivateInclude>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/cmsis-dsp/cmsis-core/Include>)
    target_compile_definitions(dsptap INTERFACE DSPTAP_FFT_CMSIS)
    message(STATUS "DspTap: float32 FFT backend = CMSIS-DSP Helium")
endif()

# ------------------------------------------------------------------------------
# Apple vDSP (Accelerate) FFT backend. When ON, the float32 real FFT in
# include/dsptap/fft.h routes through Apple's vDSP instead of Ooura — measured
# ~3x faster per transform than clang-autovectorized Ooura on Apple Silicon.
# Double is untouched (always Ooura), so the double-precision golden model is
# unaffected; on the float32 path vDSP agrees with Ooura to single-precision
# rounding (not bit-identical). Default ON for Apple, OFF everywhere else; it is
# mutually exclusive with the CMSIS backend. Force Ooura with
# -DDSPTAP_FFT_ACCELERATE=OFF.
# ------------------------------------------------------------------------------
if(APPLE AND NOT DSPTAP_FFT_CMSIS)
    set(_dsptap_fft_accel_default ON)
else()
    set(_dsptap_fft_accel_default OFF)
endif()
option(DSPTAP_FFT_ACCELERATE "Route the float32 FFT through Apple's vDSP (Accelerate) backend (Apple only)" ${_dsptap_fft_accel_default})
if(DSPTAP_FFT_ACCELERATE)
    if(NOT APPLE)
        message(FATAL_ERROR "DSPTAP_FFT_ACCELERATE requires Apple's Accelerate framework (macOS).")
    endif()
    if(DSPTAP_FFT_CMSIS)
        message(FATAL_ERROR "DSPTAP_FFT_ACCELERATE and DSPTAP_FFT_CMSIS are mutually exclusive.")
    endif()
    target_compile_definitions(dsptap INTERFACE DSPTAP_FFT_ACCELERATE)
    target_link_libraries(dsptap INTERFACE "-framework Accelerate")
    message(STATUS "DspTap: float32 FFT backend = Apple vDSP (Accelerate)")
endif()

# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------

if(PROJECT_IS_TOP_LEVEL)
    option(DSPTAP_BUILD_TESTS "Build DspTap tests" ON)
else()
    option(DSPTAP_BUILD_TESTS "Build DspTap tests" OFF)
endif()

# Warning flags for this project's own tests; never exported to consumers of
# the INTERFACE library, and never applied to third_party C.
option(DSPTAP_WERROR "Treat warnings as errors in DspTap's own targets" OFF)
add_library(dsptap_warnings INTERFACE)
target_compile_options(dsptap_warnings INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic -Wconversion -Wshadow>
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive->)
if(DSPTAP_WERROR)
    target_compile_options(dsptap_warnings INTERFACE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Werror>
        $<$<CXX_COMPILER_ID:MSVC>:/WX>)
endif()

if(DSPTAP_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
