// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Test support: the battery playing the client. Set flush-to-zero on the
// calling thread the way a client would (x86: MXCSR FTZ and DAZ; AArch64:
// FPCR.FZ; Cortex-M with an FPU: FPSCR.FZ), and restore exactly the previous
// control word on destruction, so the rest of the battery on the same thread
// is untouched. The library never touches the control register; only its
// clients (and these tests, standing in for one) do.

#pragma once

#include <cstdint>

#if defined(__x86_64__) || defined(_M_X64) || defined(__SSE_MATH__) || defined(__SSE2_MATH__)                          \
    || (defined(_M_IX86_FP) && _M_IX86_FP >= 1)
#include <xmmintrin.h>
#define TAP_DSP_FTZ_X86 1
#elif defined(_M_ARM64) && defined(_MSC_VER)
#include <intrin.h>
#define TAP_DSP_FTZ_A64_MSVC 1
#endif

namespace tap_test {

    /// The calling thread's floating-point control word, as an integer.
    inline std::uint64_t fp_control_word() noexcept {
#if defined(TAP_DSP_FTZ_X86)
        return _mm_getcsr();
#elif defined(TAP_DSP_FTZ_A64_MSVC)
        return static_cast<std::uint64_t>(_ReadStatusReg(ARM64_SYSREG(3, 3, 4, 4, 0)));
#elif defined(__aarch64__)
        std::uint64_t fpcr = 0;
        __asm__ volatile("mrs %0, fpcr" : "=r"(fpcr));
        return fpcr;
#elif defined(__ARM_FP)
        return __builtin_arm_get_fpscr();
#else
        return 0;
#endif
    }

    /// The control bits of the word only, without the sticky status flags
    /// arithmetic sets (x86 MXCSR bits 0-5; AArch32 FPSCR's NZCV, QC and
    /// cumulative flags; AArch64 keeps its flags in FPSR, not FPCR).
    inline std::uint64_t fp_control_bits() noexcept {
#if defined(TAP_DSP_FTZ_X86)
        return fp_control_word() & ~std::uint64_t{0x3F};
#elif defined(TAP_DSP_FTZ_A64_MSVC) || defined(__aarch64__)
        return fp_control_word();
#elif defined(__ARM_FP)
        return fp_control_word() & 0x07FF9F00u;
#else
        return 0;
#endif
    }

    inline void set_fp_control_word(std::uint64_t w) noexcept {
#if defined(TAP_DSP_FTZ_X86)
        _mm_setcsr(static_cast<unsigned int>(w));
#elif defined(TAP_DSP_FTZ_A64_MSVC)
        _WriteStatusReg(ARM64_SYSREG(3, 3, 4, 4, 0), static_cast<__int64>(w));
#elif defined(__aarch64__)
        __asm__ volatile("msr fpcr, %0" : : "r"(w));
#elif defined(__ARM_FP)
        __builtin_arm_set_fpscr(static_cast<unsigned int>(w));
#else
        (void)w;
#endif
    }

    /// The bits a host sets: FTZ and DAZ on x86, FZ elsewhere; 0 where none exists.
    inline constexpr std::uint64_t k_flush_bits =
#if defined(TAP_DSP_FTZ_X86)
        0x8040u;
#elif defined(TAP_DSP_FTZ_A64_MSVC) || defined(__aarch64__) || defined(__ARM_FP)
        1u << 24;
#else
        0u;
#endif

    /// Whether this target has a flush bit governing S (the soft-float M4 and
    /// double on the M4F / M33 have none; Hexagon is not modelled).
    template <typename S>
    constexpr bool has_flush_bit() noexcept {
#if defined(TAP_DSP_FTZ_X86) || defined(TAP_DSP_FTZ_A64_MSVC) || defined(__aarch64__)
        return true;
#elif defined(__ARM_FP)
        return sizeof(S) == 4 ? (__ARM_FP & 0x4) != 0 : (__ARM_FP & 0x8) != 0;
#else
        return false;
#endif
    }

    /// Sets (on = true) or clears the flush bits for its lifetime; restores the word after.
    class scoped_flush_to_zero {
      public:
        explicit scoped_flush_to_zero(bool on = true) noexcept
            : m_saved(fp_control_word()) {
            set_fp_control_word(on ? (m_saved | k_flush_bits) : (m_saved & ~k_flush_bits));
        }
        ~scoped_flush_to_zero() { set_fp_control_word(m_saved); }
        scoped_flush_to_zero(const scoped_flush_to_zero&)            = delete;
        scoped_flush_to_zero& operator=(const scoped_flush_to_zero&) = delete;

      private:
        std::uint64_t m_saved;
    };

} // namespace tap_test
