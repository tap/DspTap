// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// THE SPLIT-RADIX FINGERPRINTS (Decisions D6 and D10 of
// docs/audit-fft-and-code-smells.md): the bit patterns of
// tap::dsp::detail::split_radix_rdft's double and float outputs, forward and
// inverse, at N = 4, 16, 256, 4096 and 65536, pinned as FNV-1a-64 folds. The
// procedure, the materials and the pins live in
// tests/support/split_radix_fingerprints.h.
//
// What the pins are. Each value was measured from this engine AND from the
// vendored C it replaced (fftsg.c for double, fftsg_float.c for float, both
// at -ffp-contract=off) in one binary, and the two were equal: locally at
// 6f6f77f on x86-64 glibc (GCC 13.3.0 and clang 18.1.3) and on the four
// emulated Cortex-M legs; the parity gate's ReferenceCHasThePinnedFingerprints
// (tests/test_fft_parity_ooura.cpp) computes the same values from the C on
// every CI leg while the C is in the tree. The pins therefore carry the
// Stage 2a parity gate's promise — the engine is the C, bit for bit (D10) —
// past the gate itself.
//
// Why this target owns -ffp-contract=off. Bit identity with the C held for
// the C and the port compiled to the same sequence of IEEE operations;
// g++ contracts a*b+c into an FMA across statements by default and clang
// within one, so on FMA hardware (Apple arm64, the M55, x86 at -march) a
// default-flags build computes other bits for reasons unrelated to the
// engine's statements. MSVC: nothing is passed; /fp:precise does not
// contract (tests/CMakeLists.txt).
//
// A platform the double pin table does not name fails after printing its
// values: its libm is unmeasured, and a new platform is recorded — measured
// against the C — not skipped.

#include <cstddef>
#include <cstdint>
#include <cstdio>

#include <gtest/gtest.h>

#include "support/split_radix_fingerprints.h"
#include "tap/dsp/fft/split_radix.h"

#ifndef TAP_DSP_TEST_MAX_FFT_N
#define TAP_DSP_TEST_MAX_FFT_N 65536
#endif

namespace {

    namespace t = tap::dsp::test;

    template <typename Sample>
    using engine = tap::dsp::detail::split_radix_rdft<Sample>;

    unsigned long hi32(std::uint64_t v) {
        return static_cast<unsigned long>(v >> 32);
    }
    unsigned long lo32(std::uint64_t v) {
        return static_cast<unsigned long>(v & 0xffffffffu);
    }

    template <typename Sample>
    void expect_pinned(const char* precision, const char* platform, const t::fingerprint_pins& pins) {
        for (std::size_t i = 0; i < t::k_fingerprint_sizes.size(); ++i) {
            const std::size_t n = t::k_fingerprint_sizes[i];
            if (n > static_cast<std::size_t>(TAP_DSP_TEST_MAX_FFT_N)) {
                continue;
            }
            const t::fingerprint_pair got = t::fingerprint_of<engine<Sample>, Sample>(n);
            // Printed before any assertion so a -V log carries every host's
            // values whether or not they match (two 32-bit halves: newlib's
            // printf has no %llx).
            std::printf("[ fingerprint ] %s n=%lu forward=%08lx%08lx inverse=%08lx%08lx (%s)\n", precision,
                        static_cast<unsigned long>(n), hi32(got.forward), lo32(got.forward), hi32(got.inverse),
                        lo32(got.inverse), platform != nullptr ? platform : "no pins for this platform");
            if (platform == nullptr) {
                ADD_FAILURE() << precision << " n=" << n
                              << ": no pins for this platform's C library; record the values printed above"
                              << " (measured against the C, docs/fft-design.md) before relying on this build";
                continue;
            }
            EXPECT_NE(pins[i].forward, 0u) << precision << " n=" << n << ": unmeasured on " << platform;
            EXPECT_EQ(got.forward, pins[i].forward)
                << precision << " forward n=" << n << " on " << platform
                << ": the engine's output is not the pinned bit pattern (the vendored C's)";
            EXPECT_EQ(got.inverse, pins[i].inverse)
                << precision << " inverse n=" << n << " on " << platform
                << ": the engine's output is not the pinned bit pattern (the vendored C's)";
        }
    }

    TEST(split_radix_fingerprint, DoubleOutputIsTheVendoredCBitPattern) {
        expect_pinned<double>("double", t::k_double_pins_platform, t::k_double_pins);
    }

    TEST(split_radix_fingerprint, FloatOutputIsTheVendoredCBitPattern) {
        expect_pinned<float>("float", "every platform", t::k_float_pins);
    }

} // namespace
