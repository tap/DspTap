# Cross-compilation toolchain for Arm Cortex-M33 (bare metal, newlib +
# semihosting), executed on QEMU's MPS2+ AN505 board model. This is the
# Raspberry Pi Pico 2 W (RP2350) class of core: single-precision FPU only,
# no FP64, no MVE/Helium. The float32 profile is the profile here; anything
# double is soft-float and the on-target selection keeps it to the small
# contract tests. Ported from MuTap's cmake/arm-cortex-m33-mps2.cmake (which
# ported it from RatioTap's, which ported it from SampleRateTap's).
#
# Usage:
#   cmake -B build-m33 -DCMAKE_TOOLCHAIN_FILE=cmake/arm-cortex-m33-mps2.cmake \
#         -DCMAKE_BUILD_TYPE=MinSizeRel
# with arm-none-eabi-g++ and qemu-system-arm on PATH.
set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR arm)

set(CMAKE_C_COMPILER arm-none-eabi-gcc)
set(CMAKE_CXX_COMPILER arm-none-eabi-g++)
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

set(CMAKE_C_FLAGS_INIT "-mcpu=cortex-m33 -mthumb -mfloat-abi=hard -ffunction-sections -fdata-sections")
# -Wno-psabi: GCC otherwise emits ~26 informational "parameter passing for
# argument of type 'std::span<...>' changed in GCC 7.1" notes per leg, which
# -Werror cannot catch and which bury a real warning; the ABI note is moot in
# a single-toolchain static image.
set(CMAKE_CXX_FLAGS_INIT "${CMAKE_C_FLAGS_INIT} -Wno-psabi")

get_filename_component(_tap_dsp_platform "${CMAKE_CURRENT_LIST_DIR}/../platform" ABSOLUTE)
# Same startup as the M55 leg (Armv8-M, shared); the AN505 linker script
# places everything in the board's secure aliases (4 MB code, 4 MB data).
set(CMAKE_EXE_LINKER_FLAGS_INIT
    "--specs=rdimon.specs -nostartfiles -Wl,--gc-sections -T${_tap_dsp_platform}/mps2_an505.ld -x c ${_tap_dsp_platform}/cortexm_startup.c -x none")

set(CMAKE_CROSSCOMPILING_EMULATOR
    "qemu-system-arm;-M;mps2-an505;-nographic;-semihosting;-kernel")

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

# No Helium on the M33: the root CMakeLists defaults the CMSIS-DSP Helium FFT
# backend ON for any Generic+arm system, so pin the Ooura float32 path here
# (a plain `set` of the cache entry, which the option() then respects; the
# toolchain is processed inside project(), before that option() runs).
set(TAP_DSP_FFT_CMSIS OFF CACHE BOOL "No MVE on the Cortex-M33: Ooura float32 FFT")

# Largest transform the Stage 2a Ooura parity suite runs on this leg
# (tests/CMakeLists.txt reads it as a cache default): the plan's "parity at
# N <= 4096" on emulated targets, and what fits the data region (2^20 needs
# five 8 MB buffers). A plain cache set, so -D on the command line still wins.
set(TAP_DSP_PARITY_MAX_N 4096 CACHE STRING "Largest FFT size the Ooura parity suite runs on this leg")

# One-shot CTest mode (no argv on bare metal; see tests/CMakeLists.txt).
set(TAP_DSP_BARE_METAL ON)
