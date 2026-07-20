# Cross-compilation toolchain for Arm Cortex-M55, bare metal (newlib +
# semihosting), executed on QEMU's MPS3 AN547 board model. Ported from
# SampleRateTap's cmake/arm-cortex-m55-mps3.cmake.
#
# Usage:
#   cmake -B build-m55 -DCMAKE_TOOLCHAIN_FILE=cmake/arm-cortex-m55-mps3.cmake \
#         -DCMAKE_BUILD_TYPE=MinSizeRel
# with arm-none-eabi-g++ and qemu-system-arm on PATH.
#
# Notes:
#  - Bare metal: no std::thread (the test build adapts; see
#    tests/CMakeLists.txt and MUTAP_BARE_METAL below).
#  - The M55 FPU has no double precision, so the library's double
#    instantiations run soft-float here: correctness coverage of the
#    float32 embedded profile plus the float-tracks-double oracle check,
#    not a performance measurement.
set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR arm)

set(CMAKE_C_COMPILER arm-none-eabi-gcc)
set(CMAKE_CXX_COMPILER arm-none-eabi-g++)
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

set(CMAKE_C_FLAGS_INIT "-mcpu=cortex-m55 -mthumb -mfloat-abi=hard -ffunction-sections -fdata-sections")
set(CMAKE_CXX_FLAGS_INIT "${CMAKE_C_FLAGS_INIT}")

get_filename_component(_mutap_platform "${CMAKE_CURRENT_LIST_DIR}/../platform" ABSOLUTE)
# The startup .c is handed to the link line directly; the gcc driver
# compiles it with the same -mcpu/-mfloat-abi flags as everything else.
# `-x c` forces C compilation even under the g++ driver (which would treat
# a .c link input as C++): C guarantees the vector table's address-constant
# initializers are link-time constants, never dynamic initialization.
set(CMAKE_EXE_LINKER_FLAGS_INIT
    "--specs=rdimon.specs -nostartfiles -Wl,--gc-sections -T${_mutap_platform}/mps3_an547.ld -x c ${_mutap_platform}/armv8m_startup.c -x none")

set(CMAKE_CROSSCOMPILING_EMULATOR
    "qemu-system-arm;-M;mps3-an547;-nographic;-semihosting;-kernel")

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

# Switches the test harness to one-shot mode: a single registered CTest test
# running the whole (emulation-sized) suite, judged by gtest's summary text
# rather than the exit code, which semihosting does not reliably propagate.
set(MUTAP_BARE_METAL ON)
