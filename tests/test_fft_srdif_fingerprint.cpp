// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// THE OUTPUT FINGERPRINTS of the floating real-FFT engine, detail::srdif_rdft
// (include/tap/dsp/fft/srdif.h). For each profile and every power of two
// N = 4 … 65536 (capped by TAP_DSP_TEST_MAX_FFT_N: 4096 on the QEMU legs), one
// FNV-1a-64 over the bit pattern of every output word of a forward transform
// of a fixed input, followed by every output word of the inverse transform
// of that spectrum. The pins are in tests/support/srdif_fingerprints.h.
//
// What a pin means. The engine's arithmetic is IEEE 754 addition,
// subtraction and multiplication in a fixed order, and its tables come from
// integer arithmetic (srdif_trig: no libm anywhere in the engine), so for
// a fixed input the output bits are a function of the source alone once
// fp-contraction is excluded, which this target does (-ffp-contract=off on
// GCC / Clang / AppleClang / IntelLLVM, IntelLLVM also at -fp-model=precise;
// MSVC does not contract by default, tests/CMakeLists.txt), and every
// operation rounds to its own format (FLT_EVAL_METHOD == 0, asserted below:
// an x87 build, i386 or -mfpmath=387, evaluates in extended precision). The
// expectation is therefore ONE row for every compiler and target CI runs —
// x86-64 glibc under g++ and clang++ (either libm dispatch), Windows x64
// MSVC, macOS arm64 AppleClang, and the four Cortex-M legs under QEMU
// (arm-none-eabi-gcc: soft-float libgcc, fpv4-sp, fpv5-sp, M55) — and a
// second row would be a finding to explain, not a toolchain to record. The
// scope is those builds: a compiler this target does not name whose
// default is not IEEE-strict, MSVC on ARM64 (no leg; its default
// contraction is unverified here), or a consumer's -ffast-math or
// flush-to-zero mode can move bits the pinned input never exercises. A pin
// that moves is a numeric change to the engine: the docstring of
// fft/srdif.h and docs/fft-design.md ("The floating engine (srdif)") say
// what the numbers are and how they were measured.
//
// The engine is called directly, not through basic_real_fft (whose routing
// to its engine, byte for byte, is test_fft_routing.cpp's promise), and it
// is named explicitly, so the macOS and M55 legs, whose float default is an
// accelerated engine, pin the srdif engine too.
//
// The values are printed before anything is asserted, so a failing leg's
// log carries every value it computed; ctest swallows the output of a
// passing test, so CI runs this target verbosely as its own step
// (.github/workflows/ci.yml).

#include <cfloat>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "support/signals.h"
#include "support/srdif_fingerprints.h"
#include "tap/dsp/fft/srdif.h"

#ifndef TAP_DSP_TEST_MAX_FFT_N
#define TAP_DSP_TEST_MAX_FFT_N (1 << 20)
#endif

// The pins assume every float and double operation rounds to its own format.
#if defined(FLT_EVAL_METHOD)
static_assert(FLT_EVAL_METHOD == 0, "srdif fingerprints: evaluation in wider precision (x87) moves the output bits");
#elif !(defined(_MSC_VER) && (defined(_M_X64) || defined(_M_ARM64)))
#error "srdif fingerprints: FLT_EVAL_METHOD is not defined; confirm this target rounds each operation to its format"
#endif

namespace {

    namespace pins = tap::dsp::test::srdif_fingerprints;

    constexpr std::uint64_t k_fnv_offset = 0xcbf29ce484222325ull;
    constexpr std::uint64_t k_fnv_prime  = 0x100000001b3ull;

    /// FNV-1a-64 over each sample's bit pattern, one 32- or 64-bit word per
    /// step (the bench's fold, bench/bench_common.h).
    template <typename Sample>
    std::uint64_t fold(std::uint64_t h, const std::vector<Sample>& v) {
        using bits = std::conditional_t<sizeof(Sample) == 4, std::uint32_t, std::uint64_t>;
        for (const Sample s : v) {
            bits b = 0;
            std::memcpy(&b, &s, sizeof b);
            h ^= b;
            h *= k_fnv_prime;
        }
        return h;
    }

    /// The fixed input: 24-bit integers from xorshift32, scaled by 2^-23, so
    /// every value is exact in both profiles and the float and double inputs
    /// are the same numbers.
    template <typename Sample>
    std::vector<Sample> fingerprint_input(std::size_t n) {
        tap::dsp::test::xorshift32 rng(pins::k_seed);
        std::vector<Sample>        x(n);
        for (Sample& v : x) {
            const auto word = static_cast<std::int32_t>(rng.next_u32());
            v               = static_cast<Sample>(word >> 8) * static_cast<Sample>(0x1p-23);
        }
        return x;
    }

    template <typename Sample>
    std::uint64_t fingerprint(std::size_t n) {
        tap::dsp::detail::srdif_rdft<Sample> engine(n);
        std::vector<Sample>                  a = fingerprint_input<Sample>(n);
        engine.forward_inplace(a.data());
        std::uint64_t h = fold(k_fnv_offset, a);
        engine.inverse_inplace(a.data());
        return fold(h, a);
    }

    /// Computes and prints every size's value, then requires a row of
    /// srdif_fingerprints.h to match all of them, and names it.
    template <typename Sample>
    void expect_a_pinned_row(const char* profile) {
        constexpr std::size_t      k_cap = static_cast<std::size_t>(TAP_DSP_TEST_MAX_FFT_N);
        std::vector<std::uint64_t> got;
        for (std::size_t i = 0; i < pins::k_sizes; ++i) {
            const std::size_t n = std::size_t{4} << i;
            if (n > k_cap) {
                break;
            }
            got.push_back(fingerprint<Sample>(n));
            // %lu with casts: newlib's printf on the QEMU legs parses neither
            // %zu nor PRIx64 without __STDC_FORMAT_MACROS.
            std::printf("[ fingerprint ] srdif %s N=%lu 0x%08lx%08lx\n", profile, static_cast<unsigned long>(n),
                        static_cast<unsigned long>(got.back() >> 32),
                        static_cast<unsigned long>(got.back() & 0xffffffffu));
        }
        ASSERT_GE(got.size(), 11u) << "the sweep must reach N = 4096 on every leg";
        const pins::row* matched = nullptr;
        for (const pins::row& r : pins::k_rows) {
            const std::uint64_t* expected = std::is_same_v<Sample, float> ? r.float_rows : r.double_rows;
            bool                 all      = true;
            for (std::size_t i = 0; i < got.size(); ++i) {
                all = all && got[i] == expected[i];
            }
            if (all) {
                matched = &r;
                break;
            }
        }
        if (matched != nullptr) {
            std::printf("[ fingerprint ] srdif %s matched row: %s\n", profile, matched->name);
        }
        EXPECT_NE(matched, nullptr) << "the srdif " << profile
                                    << " output matches no pinned row: a numeric change to the engine, or a host "
                                       "that computes it differently (see the file comment)";
    }

    TEST(fft_srdif_fingerprint, FloatOutputBitsAreThePinnedRow) {
        expect_a_pinned_row<float>("float");
    }

    TEST(fft_srdif_fingerprint, DoubleOutputBitsAreThePinnedRow) {
        expect_a_pinned_row<double>("double");
    }

} // namespace
