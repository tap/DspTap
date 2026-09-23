/// @file split_radix_fingerprints.h
/// @brief The pinned output fingerprints of the split-radix FFT engine (Decision D6/D10).
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// One procedure and one table of pins, shared by the fingerprint test
// (tests/test_fft_split_radix_fingerprint.cpp) and, while the reference C is
// still in the tree, by the cross-check in tests/test_fft_parity_ooura.cpp
// that computes the same fingerprints from the C.
//
// The procedure. For each N in k_fingerprint_sizes and each precision: four
// materials in a fixed order — xorshift32 broadband noise (seed 0x9E3779B9,
// random_signal, libm-free), a unit impulse, DC (all ones), and a full-scale
// +1/-1 alternation — each taken through forward_inplace and then, on that
// forward output, inverse_inplace (unnormalized). Two FNV-1a-64 folds run
// across all four materials: one over every forward output, one over every
// inverse output, each sample folded as its IEEE bit pattern (uint64 for
// double, uint32 for float), least significant byte first. One differing bit
// anywhere changes the fold, and +0 / -0 and NaN payloads count as the bits
// they are. The on-bin tone of the old parity sweep is left out on purpose:
// tone() synthesizes through libm, so its input bits would themselves vary
// by host.
//
// Why the double pins are per platform. The engine's tables are built from
// the C library's cos / sin / atan exactly as fftsg.c builds them (Decision
// D10; split_radix.h rule 2), so a last-bit difference between two libms
// moves the double outputs at every N that reads a non-trivial twiddle. It
// moved the C's outputs identically: at 6f6f77f, the last main that carried
// the reference C, the C computed the same five distinct double sets as the
// port on every platform measured (x86-64 glibc; newlib with a
// single-precision or no FPU, i.e. the M4, M4F and M33 legs; newlib on the
// M55's double-precision FPU; x64 MSVC with the UCRT; arm64 macOS). N = 4 and 16 agree everywhere (their twiddles
// are exact). The float instantiation rounds each double libm result to
// float, and on every platform measured that absorbed the libm differences:
// one float table serves all.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <vector>

#include "signals.h"

namespace tap::dsp::test {

    /// The sizes fingerprinted: n = 4 (no real post-pass), 16, 256, 4096 (the
    /// largest the emulated legs run) and 65536 (the top of the old parity
    /// sweep). Sizes above TAP_DSP_TEST_MAX_FFT_N are not run.
    inline constexpr std::array<std::size_t, 5> k_fingerprint_sizes{4, 16, 256, 4096, 65536};

    /// FNV-1a-64 over the IEEE bit patterns of a block, LSB first.
    template <typename Sample>
    void fnv1a64_fold_bits(std::uint64_t& h, const std::vector<Sample>& block) {
        using bits_t = std::conditional_t<sizeof(Sample) == 4, std::uint32_t, std::uint64_t>;
        static_assert(sizeof(bits_t) == sizeof(Sample));
        for (const Sample v : block) {
            bits_t u = 0;
            std::memcpy(&u, &v, sizeof u);
            for (std::size_t b = 0; b < sizeof u; ++b) {
                h ^= static_cast<std::uint64_t>((u >> (8 * b)) & 0xffu);
                h *= 0x100000001b3ull;
            }
        }
    }

    /// The four materials, in fold order.
    template <typename Sample>
    std::vector<std::vector<Sample>> fingerprint_materials(std::size_t n) {
        std::vector<std::vector<Sample>> m;
        m.push_back(random_signal<Sample>(n, 0x9E3779B9u)); // broadband
        std::vector<Sample> impulse(n, Sample(0));
        impulse[0] = Sample(1);
        m.push_back(impulse);
        m.push_back(std::vector<Sample>(n, Sample(1))); // DC
        std::vector<Sample> alternating(n);
        for (std::size_t i = 0; i < n; ++i) {
            alternating[i] = (i % 2 == 0) ? Sample(1) : Sample(-1);
        }
        m.push_back(alternating);
        return m;
    }

    struct fingerprint_pair {
        std::uint64_t forward = 0;
        std::uint64_t inverse = 0;
    };

    /// The procedure above, for any engine with the constructor-from-size /
    /// forward_inplace / inverse_inplace surface (split_radix_rdft; the C's
    /// rdft behind the same surface in the parity TU).
    template <typename Engine, typename Sample>
    fingerprint_pair fingerprint_of(std::size_t n) {
        Engine           engine(n);
        fingerprint_pair f{0xcbf29ce484222325ull, 0xcbf29ce484222325ull};
        for (std::vector<Sample> x : fingerprint_materials<Sample>(n)) {
            engine.forward_inplace(x.data());
            fnv1a64_fold_bits(f.forward, x);
            engine.inverse_inplace(x.data());
            fnv1a64_fold_bits(f.inverse, x);
        }
        return f;
    }

    /// One platform's pins, indexed like k_fingerprint_sizes. A zero pair is
    /// a size never measured there (the emulated legs stop at 4096).
    using fingerprint_pins = std::array<fingerprint_pair, k_fingerprint_sizes.size()>;

    /// float: every platform measured (see the file comment).
    inline constexpr fingerprint_pins k_float_pins{{
        {0xef9b60ebf0e8f587ull, 0x7d9b9e4480f66924ull}, // 4
        {0xa222e5c5b80bba3bull, 0x3a182210d01399a7ull}, // 16
        {0x1dc8ccaa8f8db08dull, 0x50306aa83541329dull}, // 256
        {0x6723317886e7c1f2ull, 0x7a318e19506b1c24ull}, // 4096
        {0x141744d674ac1c64ull, 0x8c9d260bedee7129ull}, // 65536
    }};

    /// double, per C library. The selection is by the macros that identify
    /// the libm the engine's tables are built from; a platform none of them
    /// names has no pins (k_double_pins_platform is nullptr).
#if defined(__GLIBC__) && defined(__x86_64__)
    inline constexpr const char*      k_double_pins_platform = "x86-64 glibc";
    inline constexpr fingerprint_pins k_double_pins{{
        {0xf014008109d04382ull, 0xa4e32575aed1c801ull}, // 4
        {0xea32e166ad289557ull, 0x4a89b6e54a2b629eull}, // 16
        {0x4e70f1158efcfe3dull, 0xb5d51404f905fea6ull}, // 256
        {0x753d7cfca461e82bull, 0x8c9e917a44067ff2ull}, // 4096
        {0xb7909fad30a275afull, 0xe72fb50d719ef181ull}, // 65536
    }};
#elif defined(__NEWLIB__) && defined(__arm__) && defined(__ARM_FP) && (__ARM_FP & 0x8)
    inline constexpr const char*      k_double_pins_platform = "Arm newlib, double-precision FPU (M55)";
    inline constexpr fingerprint_pins k_double_pins{{
        {0xf014008109d04382ull, 0xa4e32575aed1c801ull}, // 4
        {0xea32e166ad289557ull, 0x4a89b6e54a2b629eull}, // 16
        {0x5c6c45c5898a558aull, 0x43318e1983699fd3ull}, // 256
        {0xc42e9fb613c1c2c9ull, 0x81d554f525c82f66ull}, // 4096
        {0, 0},                                         // 65536: not run
    }};
#elif defined(__NEWLIB__) && defined(__arm__)
    inline constexpr const char*      k_double_pins_platform = "Arm newlib, software double (M4, M4F, M33)";
    inline constexpr fingerprint_pins k_double_pins{{
        {0xf014008109d04382ull, 0xa4e32575aed1c801ull}, // 4
        {0xea32e166ad289557ull, 0x4a89b6e54a2b629eull}, // 16
        {0xe59fe4d3e6326eb3ull, 0x43318e1983699fd3ull}, // 256
        {0x58c61a3b29d4b55eull, 0x2d0482b12bc4e246ull}, // 4096
        {0, 0},                                         // 65536: not run
    }};
#elif defined(_MSC_VER) && defined(_M_X64)
    inline constexpr const char*      k_double_pins_platform = "x64 MSVC, UCRT";
    inline constexpr fingerprint_pins k_double_pins{{
        {0xf014008109d04382ull, 0xa4e32575aed1c801ull}, // 4
        {0xea32e166ad289557ull, 0x4a89b6e54a2b629eull}, // 16
        {0xcabb356393c2dac7ull, 0x5c0eeb0579868a88ull}, // 256
        {0x4e0db6e49731d286ull, 0xbcbcc885236c3585ull}, // 4096
        {0xb02d408a57e2ebdeull, 0x288c58dad9580f4bull}, // 65536
    }};
#elif defined(__APPLE__) && defined(__aarch64__)
    inline constexpr const char*      k_double_pins_platform = "arm64 macOS";
    inline constexpr fingerprint_pins k_double_pins{{
        {0xf014008109d04382ull, 0xa4e32575aed1c801ull}, // 4
        {0xea32e166ad289557ull, 0x4a89b6e54a2b629eull}, // 16
        {0x1a682072ef54d745ull, 0xa9e60433894b7013ull}, // 256
        {0x7909ca3d02d78252ull, 0x4a9a5c439158323dull}, // 4096
        {0xbae4e82f8d79d549ull, 0xbfa3a82a87bbe16dull}, // 65536
    }};
#else
    inline constexpr const char*      k_double_pins_platform = nullptr;
    inline constexpr fingerprint_pins k_double_pins{};
#endif

} // namespace tap::dsp::test
