# Cross-compilation toolchain for Arm Cortex-M4 (bare metal, newlib +
# semihosting), executed on QEMU's MPS2 AN386 board model. Two flavours,
# selected by the TAP_DSP_M4_FPU cache variable:
#   OFF (default)  -mfloat-abi=soft: no FPU at all. Every float and double
#                  operation is a library call — the fixed-point profiles'
#                  reason to exist, and the leg where any accidental double
#                  in a float or fixed path shows up in the emulation budget.
#   ON             -mfpu=fpv4-sp-d16 -mfloat-abi=hard (M4F): single-precision
#                  hardware FPU, no FP64, no DSP-extension SIMD beyond what
#                  Armv7E-M carries (SMLALD etc. are present in both flavours).
#
# Usage:
#   cmake -B build-m4 -DCMAKE_TOOLCHAIN_FILE=cmake/arm-cortex-m4-mps2.cmake \
#         -DCMAKE_BUILD_TYPE=MinSizeRel [-DTAP_DSP_M4_FPU=ON]
# with arm-none-eabi-g++ and qemu-system-arm on PATH.
set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR arm)

set(CMAKE_C_COMPILER arm-none-eabi-gcc)
set(CMAKE_CXX_COMPILER arm-none-eabi-g++)
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

set(TAP_DSP_M4_FPU OFF CACHE BOOL "Cortex-M4F: fpv4-sp-d16 hard-float ABI (OFF = soft-float, no FPU)")
# Toolchain files are re-read inside try_compile projects, which do not see
# the parent cache; forward the switch so both reads pick the same flags.
set(CMAKE_TRY_COMPILE_PLATFORM_VARIABLES TAP_DSP_M4_FPU)

if(TAP_DSP_M4_FPU)
    set(_tap_dsp_m4_float "-mfpu=fpv4-sp-d16 -mfloat-abi=hard")
else()
    set(_tap_dsp_m4_float "-mfloat-abi=soft")
endif()
set(CMAKE_C_FLAGS_INIT "-mcpu=cortex-m4 -mthumb ${_tap_dsp_m4_float} -ffunction-sections -fdata-sections")
# -Wno-psabi: GCC otherwise emits ~26 informational "parameter passing for
# argument of type 'std::span<...>' changed in GCC 7.1" notes per leg, which
# -Werror cannot catch and which bury a real warning; the ABI note is moot in
# a single-toolchain static image.
set(CMAKE_CXX_FLAGS_INIT "${CMAKE_C_FLAGS_INIT} -Wno-psabi")

get_filename_component(_tap_dsp_platform "${CMAKE_CURRENT_LIST_DIR}/../platform" ABSOLUTE)
# Shared Cortex-M startup (the MSPLIM write is compiled out on Armv7E-M);
# the AN386 linker script places code in ZBT SSRAM1 and data in SSRAM2/3
# (4 MB each). `-x c` forces C compilation of the startup under the g++
# driver: C guarantees the vector table's address-constant initializers are
# link-time constants, never dynamic initialization.
set(CMAKE_EXE_LINKER_FLAGS_INIT
    "--specs=rdimon.specs -nostartfiles -Wl,--gc-sections -T${_tap_dsp_platform}/mps2_an386.ld -x c ${_tap_dsp_platform}/cortexm_startup.c -x none")

set(CMAKE_CROSSCOMPILING_EMULATOR
    "qemu-system-arm;-M;mps2-an386;-nographic;-semihosting;-kernel")

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

# No Helium on the M4: the compiler defines no __ARM_FEATURE_MVE for this
# -mcpu, so the root CMakeLists.txt's MVE-F check leaves TAP_DSP_FFT_CMSIS OFF
# (the split-radix float32 engine) without a pin here. Until that check the
# default was ON for every Generic+arm system and this file pinned it OFF.

# Largest transform the FFT test sweeps run on this leg (tests/CMakeLists.txt
# reads it as a cache default): what fits the data region (2^20 needs five
# 8 MB buffers). Introduced as TAP_DSP_PARITY_MAX_N for the Stage 2a parity
# sweep, renamed when that gate was retired (Decision D6). A plain cache set,
# so -D on the command line still wins.
set(TAP_DSP_TEST_MAX_FFT_N 4096 CACHE STRING "Largest FFT size the test sweeps run on this leg")

# One-shot CTest mode (no argv on bare metal; see tests/CMakeLists.txt).
set(TAP_DSP_BARE_METAL ON)
