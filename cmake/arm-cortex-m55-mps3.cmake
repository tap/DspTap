# Cross-compilation toolchain for Arm Cortex-M55, bare metal (newlib +
# semihosting), executed on QEMU's MPS3 AN547 board model. Ported from
# MuTap's cmake/arm-cortex-m55-mps3.cmake (itself from SampleRateTap's).
# The only leg with MVE/Helium: -mcpu=cortex-m55 -mfloat-abi=hard defines
# __ARM_FEATURE_MVE = 3 (integer + floating-point MVE), so the root
# CMakeLists.txt's MVE-F check defaults the CMSIS-DSP float32 FFT backend ON
# here, and its parity suite runs under emulation.
#
# Usage:
#   cmake -B build-m55 -DCMAKE_TOOLCHAIN_FILE=cmake/arm-cortex-m55-mps3.cmake \
#         -DCMAKE_BUILD_TYPE=MinSizeRel
# with arm-none-eabi-g++ and qemu-system-arm on PATH.
#
# Notes:
#  - Bare metal: no std::thread (the test build adapts; see
#    tests/CMakeLists.txt and TAP_DSP_BARE_METAL below).
#  - The double instantiations are correctness coverage of the golden model
#    on a 32-bit target plus the float-tracks-double oracle checks, not a
#    performance measurement.
set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR arm)

set(CMAKE_C_COMPILER arm-none-eabi-gcc)
set(CMAKE_CXX_COMPILER arm-none-eabi-g++)
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

set(CMAKE_C_FLAGS_INIT "-mcpu=cortex-m55 -mthumb -mfloat-abi=hard -ffunction-sections -fdata-sections")
# -Wno-psabi: GCC otherwise emits ~26 informational "parameter passing for
# argument of type 'std::span<...>' changed in GCC 7.1" notes per leg, which
# -Werror cannot catch and which bury a real warning; the ABI note is moot in
# a single-toolchain static image.
set(CMAKE_CXX_FLAGS_INIT "${CMAKE_C_FLAGS_INIT} -Wno-psabi")

get_filename_component(_tap_dsp_platform "${CMAKE_CURRENT_LIST_DIR}/../platform" ABSOLUTE)
# The startup .c is handed to the link line directly; the gcc driver
# compiles it with the same -mcpu/-mfloat-abi flags as everything else.
# `-x c` forces C compilation even under the g++ driver (which would treat
# a .c link input as C++): C guarantees the vector table's address-constant
# initializers are link-time constants, never dynamic initialization.
set(CMAKE_EXE_LINKER_FLAGS_INIT
    "--specs=rdimon.specs -nostartfiles -Wl,--gc-sections -T${_tap_dsp_platform}/mps3_an547.ld -x c ${_tap_dsp_platform}/cortexm_startup.c -x none")

set(CMAKE_CROSSCOMPILING_EMULATOR
    "qemu-system-arm;-M;mps3-an547;-nographic;-semihosting;-kernel")

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

# Largest transform the FFT test sweeps run on this leg (tests/CMakeLists.txt
# reads it as a cache default): what fits the data region (2^20 needs five
# 8 MB buffers). Introduced as TAP_DSP_PARITY_MAX_N for the Stage 2a parity
# sweep, renamed when that gate was retired (Decision D6). A plain cache set,
# so -D on the command line still wins.
set(TAP_DSP_TEST_MAX_FFT_N 4096 CACHE STRING "Largest FFT size the test sweeps run on this leg")

# Switches the test harness to one-shot mode: a single registered CTest test
# running the whole (emulation-sized) suite, judged by gtest's summary text
# rather than the exit code, which semihosting does not reliably propagate.
set(TAP_DSP_BARE_METAL ON)
