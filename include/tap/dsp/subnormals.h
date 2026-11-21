/// @file subnormals.h
/// @brief Whether the calling thread flushes subnormal results to zero, per sample type and target.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// A recursive (IIR) stage's state decays into the subnormal range once its
// input goes silent. Two things then go wrong: subnormal arithmetic is slow on
// x86 and on some A-class Arm cores (the classic CPU spike on silence), and a
// state can stick at a small multiple of the smallest subnormal instead of
// reaching zero. FIR stages have neither problem: their delay lines empty
// after one window.
//
// Flushing subnormals is the HOST's job (SampleRateTap's IIR plan, I9,
// decided by the maintainer): the control bit belongs to the thread, most
// audio hosts already set it on their audio threads, and a library that set
// it would have to restore it on every call. No DspTap stage reads-and-writes
// or restores the floating-point control register. This header only lets a
// stage, a host and a test ASK, so a recursive stage can state the
// requirement as a debug-build precondition:
//
//   TAP_EXPECTS(tap::dsp::subnormals_flushed<S>())   (where a bit governs S)
//
// What the host sets, per target (the plan's table, I9):
//
//   | target                                     | bit                              | governs            |
//   |--------------------------------------------|----------------------------------|--------------------|
//   | x86-64 (and x86 with SSE math)             | MXCSR FTZ (bit 15); DAZ (bit 6)  | float and double   |
//   | AArch64 (Apple silicon, Arm Linux, WoA)    | FPCR.FZ (bit 24)                 | float and double   |
//   | AArch32 with an FPU (Cortex-M4F, M33, M55) | FPSCR.FZ (bit 24); FPDSCR too,   | float; double only |
//   |                                            | so exception handlers inherit it | where the FPU has  |
//   |                                            |                                  | double (M55)       |
//   | Cortex-M without an FPU (M4), double on    | none: software float             | -                  |
//   | the M4F and M33, Hexagon                   |                                  |                    |
//
// How a host sets it, on each audio thread before it runs (x86:
// `_mm_setcsr(_mm_getcsr() | 0x8040)`; AArch64: `msr fpcr` with bit 24;
// Cortex-M: `__builtin_arm_set_fpscr(__builtin_arm_get_fpscr() | (1u << 24))`
// or CMSIS's `__set_FPSCR`, at reset). DAZ (x86) is recommended for speed and
// not required by the query: FTZ alone is what makes a decaying state reach
// exact zero.
//
// Contract, as numbers:
//   - has_subnormal_flush_control<S>() is true exactly where a row above names
//     a bit for S, decided at compile time from the target macros (x86:
//     __SSE_MATH__ / __SSE2_MATH__, or _M_X64 under MSVC; AArch64:
//     __aarch64__ / _M_ARM64; AArch32: __ARM_FP bit 2 for float, bit 3 for
//     double). Hexagon and soft-float targets are false.
//   - subnormals_flushed<S>() reads the bit (one register read, noexcept,
//     no write) where the trait is true, and returns false where it is not:
//     nothing flushes S's subnormals on such a target, so a stage's
//     honest limit applies there (the residue a recursive stage documents).
//   - Untested targets, stated: MSVC on ARM64 (no CI leg) reads FPCR through
//     _ReadStatusReg; x86 with x87 float math reports false.
//
// TAP_DSP_ALLOW_SUBNORMALS, defined by a consumer before any tap/dsp header,
// compiles the recursive stages' precondition out (the plan's placeholder
// name), for a host that deliberately runs without FTZ and accepts the
// documented residue. Like TAP_EXPECTS it expands inside inline functions,
// so it must be the same in every translation unit of an image.

#pragma once

#include <cstdint>
#include <type_traits>

#include "tap/dsp/detail/expects.h"

#if defined(__x86_64__) || defined(_M_X64) || defined(__SSE_MATH__) || defined(__SSE2_MATH__)                          \
    || (defined(_M_IX86_FP) && _M_IX86_FP >= 1)
#include <xmmintrin.h>
#define TAP_DSP_FTZ_X86 1
#elif defined(_M_ARM64) && defined(_MSC_VER)
#include <intrin.h>
#define TAP_DSP_FTZ_A64_MSVC 1
#endif

namespace tap::dsp {

    namespace detail {

        template <typename S>
        consteval bool subnormal_flush_controlled() {
            if constexpr (std::is_same_v<S, float>) {
#if defined(_M_X64) || defined(__SSE_MATH__) || (defined(_M_IX86_FP) && _M_IX86_FP >= 1)
                return true;
#elif defined(__aarch64__) || defined(_M_ARM64)
                return true;
#elif defined(__ARM_FP) && (__ARM_FP & 0x4)
                return true;
#else
                return false;
#endif
            }
            else if constexpr (std::is_same_v<S, double>) {
#if defined(_M_X64) || defined(__SSE2_MATH__) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
                return true;
#elif defined(__aarch64__) || defined(_M_ARM64)
                return true;
#elif defined(__ARM_FP) && (__ARM_FP & 0x8)
                return true;
#else
                return false;
#endif
            }
            else {
                return false; // integer formats have no subnormals
            }
        }

        /// The flush bit of the calling thread's control register, where one exists.
        inline bool read_flush_bit() noexcept {
#if defined(TAP_DSP_FTZ_X86)
            return (_mm_getcsr() & 0x8000u) != 0;
#elif defined(TAP_DSP_FTZ_A64_MSVC)
            return (_ReadStatusReg(ARM64_SYSREG(3, 3, 4, 4, 0)) & (1ull << 24)) != 0; // FPCR
#elif defined(__aarch64__)
            std::uint64_t fpcr = 0;
            __asm__ volatile("mrs %0, fpcr" : "=r"(fpcr));
            return (fpcr & (1ull << 24)) != 0;
#elif defined(__ARM_FP)
            return (__builtin_arm_get_fpscr() & (1u << 24)) != 0;
#else
            return false;
#endif
        }

    } // namespace detail

    /// True where a control bit flushes S's subnormal results on this target (file header table).
    template <typename S>
    constexpr bool has_subnormal_flush_control() noexcept {
        return detail::subnormal_flush_controlled<S>();
    }

    /// Whether the calling thread currently flushes S's subnormal results to zero.
    /// Reads, never writes; false wherever has_subnormal_flush_control<S>() is false.
    template <typename S>
    [[nodiscard]] inline bool subnormals_flushed() noexcept {
        if constexpr (has_subnormal_flush_control<S>()) {
            return detail::read_flush_bit();
        }
        else {
            return false;
        }
    }

} // namespace tap::dsp

/// The recursive stages' precondition (subnormals.h, file header): compiled
/// out by TAP_DSP_ALLOW_SUBNORMALS, and nothing where no bit governs S.
#if defined(TAP_DSP_ALLOW_SUBNORMALS)
#define TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED(S) ((void)0)
#else
#define TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED(S)                                                                          \
    do {                                                                                                               \
        if constexpr (::tap::dsp::has_subnormal_flush_control<S>()) {                                                  \
            TAP_EXPECTS(::tap::dsp::subnormals_flushed<S>());                                                          \
        }                                                                                                              \
    } while (false)
#endif
