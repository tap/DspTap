/// @file subnormals.h
/// @brief Whether subnormal results are flushed to zero on the calling thread: the host's job, a recursive kernel's
/// precondition.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// After its input goes silent, a floating-point recursive state decays into
// the subnormal range, and two things go wrong (SampleRateTap's
// docs/IIR_PLAN.md, I9): subnormal arithmetic is slow on x86 and on some
// A-class Arm cores, and on any target a first-order section with 1/2 < a <
// 1 sticks at k times the smallest subnormal, k <= 1 / (2 (1 - a)), because
// under round to nearest a times it rounds back to itself. FIR stages have
// neither problem (their delay lines empty after one window); the recursive
// kernels (iir_kernels.h) have both.
//
// Flushing subnormals is the host's job (I9, the maintainer's decision): the
// bit belongs to the thread, most audio hosts already set it on their audio
// threads, and a library that set it would have to restore it on every call.
// No kernel here reads-and-writes or restores the floating-point control
// register; each recursive kernel only READS it, as a debug-build
// precondition (TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED: TAP_EXPECTS, so an
// assert in a Debug build and nothing in Release), on targets where a bit
// governs the sample type's arithmetic, and compiles nothing elsewhere.
// Defining TAP_DSP_ALLOW_SUBNORMALS before any tap/dsp header compiles the
// check out (for a host that accepts the honest limit below, and for the
// tests that pin it); like TAP_EXPECTS it must be the same in every
// translation unit of an image (an ODR rule: it expands inside inline
// functions).
//
// What the host sets, per target (k_subnormal_control says which types a bit
// governs on the target being compiled for):
//
//   | target                         | what the host sets                        | governs           |
//   |--------------------------------|-------------------------------------------|-------------------|
//   | x86-64                         | MXCSR FTZ (bit 15) and DAZ (bit 6):       | float and double  |
//   |                                | _mm_setcsr(_mm_getcsr() | 0x8040)         | (SSE)             |
//   |                                | (<xmmintrin.h>; GCC, Clang, MSVC), on     |                   |
//   |                                | each audio thread before it runs          |                   |
//   | AArch64 (Apple silicon, Arm    | FPCR.FZ (bit 24): msr fpcr (inline asm    | float and double  |
//   | Linux, Windows on Arm)         | under GCC and Clang, _WriteStatusReg      |                   |
//   |                                | under MSVC), on each audio thread         |                   |
//   | Cortex-M with an FPU (M4F,     | FPSCR.FZ (bit 24), and FPDSCR             | float on the M4F  |
//   | M33, M55), AArch32             | (0xE000EF3C) so exception handlers, where | and M33; float    |
//   |                                | MCU audio often runs, start with it:      | and double on the |
//   |                                | __builtin_arm_set_fpscr (or CMSIS's       | M55               |
//   |                                | __set_FPSCR), at reset                    |                   |
//   | Cortex-M without an FPU (M4),  | nothing: software float, IEEE gradual     | nothing           |
//   | double on the M4F and M33      | underflow                                 |                   |
//   | Hexagon                        | nothing: USR has no flush-to-zero field   | nothing           |
//   |                                | (the vendor's V68 Programmer's Reference  |                   |
//   |                                | Manual, 80-N2040-46 Rev. B, Table 2-6;    |                   |
//   |                                | its floating point "defined to support    |                   |
//   |                                | the IEEE floating-point standard",        |                   |
//   |                                | 4.3.4), and QEMU models none (IIR-0)      |                   |
//
// Honest limits, for a host that does not set it or a target with no bit
// (pinned by tests/test_iir_kernels.cpp and tests/test_iir_subnormals.cpp):
//   - CPU: on x86 and some A-class Arm cores silence after signal runs 9 to
//     11 times slower once the states go subnormal (IIR-0, x86-64);
//   - output: the state sticks at k times the smallest subnormal, so the
//     output after a flush window is a residue of at most a few multiples of
//     the smallest subnormal (about 1e-45 in float, 1e-323 in double), not
//     silence.
//
// A consumer built with -ffast-math (or /fp:fast, or icx's default
// -fp-model=fast) starts with FTZ and DAZ already on on x86 (GCC links
// crtfastmath.o into the executable), which does what this header asks of a
// host; the rest of what -ffast-math does to these contracts is README.md's
// honest limit.

#pragma once

#include <cstdint>
#include <type_traits>

#include "tap/dsp/detail/expects.h"
#include "tap/dsp/sample_traits.h"

#if defined(__x86_64__) || defined(_M_X64)
#include <xmmintrin.h>
#define TAP_DSP_FP_CONTROL_X86_64 1
#elif defined(__aarch64__) || defined(_M_ARM64)
#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#endif
#define TAP_DSP_FP_CONTROL_AARCH64 1
#elif defined(__arm__) && defined(__ARM_FP)
#define TAP_DSP_FP_CONTROL_AARCH32 1
#endif

namespace tap::dsp {

    // ANCHOR: subnormal_control
    /// Whether a floating-point control bit of this target governs subnormal
    /// results of S's arithmetic (the file header's table): true for float
    /// and double on x86-64 and AArch64, for float on 32-bit Arm with a
    /// single-precision FPU and for double with a double-precision one; false
    /// for every fixed-point format, on Hexagon and with software float.
    template <sample_type S>
    inline constexpr bool k_subnormal_control =
#if defined(TAP_DSP_FP_CONTROL_X86_64) || defined(TAP_DSP_FP_CONTROL_AARCH64)
        std::is_floating_point_v<S>;
#elif defined(TAP_DSP_FP_CONTROL_AARCH32)
        // __ARM_FP: bit 2 single precision in hardware, bit 3 double.
        (std::is_same_v<S, float> && (__ARM_FP & 4) != 0) || (std::is_same_v<S, double> && (__ARM_FP & 8) != 0);
#else
        false;
#endif
    // ANCHOR_END: subnormal_control

    /// The control bits a host sets to flush subnormals on this target (0
    /// where there are none): MXCSR FTZ | DAZ on x86-64, FPCR.FZ on AArch64,
    /// FPSCR.FZ on 32-bit Arm.
    inline constexpr std::uint64_t k_flush_to_zero_bits =
#if defined(TAP_DSP_FP_CONTROL_X86_64)
        0x8040u;
#elif defined(TAP_DSP_FP_CONTROL_AARCH64) || defined(TAP_DSP_FP_CONTROL_AARCH32)
        std::uint64_t{1} << 24;
#else
        0u;
#endif

    /// The calling thread's floating-point control register as this target
    /// names it (MXCSR, FPCR, FPSCR), read without side effects; 0 where the
    /// target has none. Never written by anything in tap/dsp.
    inline std::uint64_t floating_point_control() noexcept {
#if defined(TAP_DSP_FP_CONTROL_X86_64)
        return _mm_getcsr();
#elif defined(TAP_DSP_FP_CONTROL_AARCH64)
#if defined(_MSC_VER) && !defined(__clang__)
        return static_cast<std::uint64_t>(_ReadStatusReg(ARM64_FPCR));
#else
        std::uint64_t v = 0;
        __asm__ volatile("mrs %0, fpcr" : "=r"(v));
        return v;
#endif
#elif defined(TAP_DSP_FP_CONTROL_AARCH32)
        return __builtin_arm_get_fpscr();
#else
        return 0;
#endif
    }

    // ANCHOR: subnormals_flushed
    /// Whether S's subnormal results are flushed to zero on the calling
    /// thread right now: the control bits read where k_subnormal_control<S>;
    /// true for a fixed-point S (it has no subnormals to flush); false for a
    /// floating S that no bit governs (software float, Hexagon), where the
    /// header's honest limit applies.
    template <sample_type S>
    inline bool subnormals_flushed() noexcept {
        if constexpr (sample_traits<S>::k_is_fixed_point) {
            return true;
        }
        else if constexpr (!k_subnormal_control<S>) {
            return false;
        }
        else {
            return (floating_point_control() & k_flush_to_zero_bits) == k_flush_to_zero_bits;
        }
    }
    // ANCHOR_END: subnormals_flushed

} // namespace tap::dsp

// ANCHOR: expects_flushed
/// The recursive kernels' debug-build precondition: TAP_EXPECTS that S's
/// subnormals are flushed, where a bit governs S's arithmetic; nothing
/// elsewhere, nothing in a release build, nothing under
/// TAP_DSP_ALLOW_SUBNORMALS.
#if defined(TAP_DSP_ALLOW_SUBNORMALS)
#define TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED(S) static_cast<void>(0)
#else
#define TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED(S)                                                                          \
    do {                                                                                                               \
        if constexpr (::tap::dsp::k_subnormal_control<S>) {                                                            \
            TAP_EXPECTS(::tap::dsp::subnormals_flushed<S>());                                                          \
        }                                                                                                              \
    } while (false)
#endif
// ANCHOR_END: expects_flushed
