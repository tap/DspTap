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
set(CMAKE_CXX_FLAGS_INIT "${CMAKE_C_FLAGS_INIT}")

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

# No Helium on the M4: the root CMakeLists defaults the CMSIS-DSP Helium FFT
# backend ON for any Generic+arm system, so pin the Ooura float32 path here
# (a plain `set` of the cache entry, which the option() then respects; the
# toolchain is processed inside project(), before that option() runs).
set(TAP_DSP_FFT_CMSIS OFF CACHE BOOL "No MVE on the Cortex-M4: Ooura float32 FFT")

# One-shot CTest mode (no argv on bare metal; see tests/CMakeLists.txt).
set(TAP_DSP_BARE_METAL ON)
