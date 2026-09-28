# SPDX-License-Identifier: MIT
# Copyright 2026 Timothy Place and the DspTap contributors.
# Adapted from MuTap's cmake/hexagon-linux-musl.cmake (MIT, MuTap contributors).
#
# Cross-compile DspTap for Qualcomm Hexagon: triple hexagon-unknown-linux-musl,
# built with the CodeLinaro "toolchain for hexagon" (clang + musl sysroot + LLVM
# runtimes). Point HEXAGON_TOOLCHAIN_ROOT (cache variable or environment) at
# the unpacked clang+llvm-*-cross-hexagon-unknown-linux-musl/x86_64-linux-gnu
# directory. Used by bench.yml's `hexagon` key (the instruction-count ratchet
# under qemu-hexagon user-mode emulation); MuTap's Hexagon CI leg and ratchet
# use the same toolchain, flags and QEMU, so the two repositories' Hexagon
# counts are comparable.
#
# A hosted Linux target, not bare metal: the bench binaries print through
# stdio and exit normally, and are linked statically (first-class with musl)
# so the emulator needs no sysroot path.
#
# Flags (MuTap's, kept identical):
#   -mv68                   the newest revision the shipped sysroot libraries
#                           and qemu-hexagon agree on.
#   -mhvx -mhvx-length=128b the 128-byte HVX unit, for auto-vectorization.

set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR Hexagon)

if(NOT DEFINED HEXAGON_TOOLCHAIN_ROOT OR HEXAGON_TOOLCHAIN_ROOT STREQUAL "")
    set(HEXAGON_TOOLCHAIN_ROOT "$ENV{HEXAGON_TOOLCHAIN_ROOT}")
endif()
if(HEXAGON_TOOLCHAIN_ROOT STREQUAL "")
    message(FATAL_ERROR
        "Set HEXAGON_TOOLCHAIN_ROOT (cache variable or environment) to the unpacked CodeLinaro "
        "hexagon toolchain (clang+llvm-*-cross-hexagon-unknown-linux-musl/x86_64-linux-gnu)")
endif()
# try_compile projects re-read this file without the cache: pass the root on.
set(ENV{HEXAGON_TOOLCHAIN_ROOT} "${HEXAGON_TOOLCHAIN_ROOT}")

set(TAP_DSP_HEXAGON_TRIPLE hexagon-unknown-linux-musl)
if(EXISTS "${HEXAGON_TOOLCHAIN_ROOT}/bin/${TAP_DSP_HEXAGON_TRIPLE}-clang++")
    set(CMAKE_C_COMPILER "${HEXAGON_TOOLCHAIN_ROOT}/bin/${TAP_DSP_HEXAGON_TRIPLE}-clang")
    set(CMAKE_CXX_COMPILER "${HEXAGON_TOOLCHAIN_ROOT}/bin/${TAP_DSP_HEXAGON_TRIPLE}-clang++")
else()
    set(CMAKE_C_COMPILER "${HEXAGON_TOOLCHAIN_ROOT}/bin/clang")
    set(CMAKE_CXX_COMPILER "${HEXAGON_TOOLCHAIN_ROOT}/bin/clang++")
    set(CMAKE_C_COMPILER_TARGET ${TAP_DSP_HEXAGON_TRIPLE})
    set(CMAKE_CXX_COMPILER_TARGET ${TAP_DSP_HEXAGON_TRIPLE})
endif()

set(TAP_DSP_HEXAGON_ARCH_FLAGS "-mv68 -mhvx -mhvx-length=128b")
set(CMAKE_C_FLAGS_INIT "${TAP_DSP_HEXAGON_ARCH_FLAGS}")
set(CMAKE_CXX_FLAGS_INIT "${TAP_DSP_HEXAGON_ARCH_FLAGS}")
# --eh-frame-hdr: not implied for -static links, but the unwinder needs
# PT_GNU_EH_FRAME to find the exception tables.
set(CMAKE_EXE_LINKER_FLAGS_INIT "-static -Wl,--eh-frame-hdr")

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

find_program(TAP_DSP_QEMU_HEXAGON NAMES qemu-hexagon-static qemu-hexagon)
if(TAP_DSP_QEMU_HEXAGON)
    set(CMAKE_CROSSCOMPILING_EMULATOR "${TAP_DSP_QEMU_HEXAGON}")
endif()
