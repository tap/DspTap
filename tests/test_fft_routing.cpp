// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// THE ROUTING PROOF for the Stage 2b flip (docs/audit-fft-and-code-smells.md,
// Part 3): tap::dsp::basic_real_fft<double> is detail::split_radix_rdft<double>
// and, where the build selected no float backend, basic_real_fft<float> is
// detail::split_radix_rdft<float> — byte for byte, forward and inverse.
//
// Why this is a test and not a reading of fft.h. The bit-identity gate
// (tests/test_fft_parity_ooura.cpp) proves the ENGINE is the vendored C; it
// deliberately does not go through basic_real_fft, because on the M55 and
// macOS legs basic_real_fft<float> is a backend, not the engine. This file
// closes the remaining gap: that the consumer-facing class adds nothing
// between the caller and the engine (no scale, no copy that could round, no
// stale table) on the profiles that route to it. Together the two files say
// what the Stage 2b PR claims: every consumer's output through
// basic_real_fft is the C's, bit for bit.
//
// memcmp, not EXPECT_EQ: == calls +0.0 and -0.0 equal and any NaN unequal to
// itself; the claim is about bits. Both sides receive identical input bits,
// and the inverse is judged on the engine's own forward output, so the two
// directions are pinned independently.
//
// The float half is compiled out under TAP_DSP_FFT_CMSIS / TAP_DSP_FFT_ACCELERATE
// (fft.h defines TAP_DSP_FFT_FLOAT_BACKEND for either): there
// basic_real_fft<float> is the backend by design and its contract is pinned
// to the engine at float epsilon by tests/test_fft_backend.cpp instead.

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "support/signals.h"
#include "tap/dsp/fft.h"
#include "tap/dsp/fft/split_radix.h"

namespace {

    // Sizes small enough for the QEMU legs (Part 10) and wide enough to cross
    // every code path the engine dispatches on: n = 4 (no real post-pass,
    // makect never called), 8 (the smallest post-pass), the certified
    // geometries 512 and 2048, and 4096 (the largest the emulated legs run).
    constexpr std::size_t k_sizes[] = {4, 8, 64, 512, 2048, 4096};

    template <typename Sample>
    std::vector<std::vector<Sample>> materials(std::size_t n) {
        std::vector<std::vector<Sample>> m;
        m.push_back(tap::dsp::test::random_signal<Sample>(n, 0x9E3779B9u)); // broadband
        m.push_back(tap::dsp::test::tone<Sample>(n, static_cast<double>(std::max<std::size_t>(1, n / 8)), 0.5, 0.3));
        std::vector<Sample> impulse(n, Sample(0));
        impulse[0] = Sample(1);
        m.push_back(impulse);
        m.push_back(std::vector<Sample>(n, Sample(1))); // dc
        return m;
    }

    template <typename Sample>
    void expect_routed_to_the_engine() {
        for (const std::size_t n : k_sizes) {
            tap::dsp::basic_real_fft<Sample>           fft(n);
            tap::dsp::detail::split_radix_rdft<Sample> engine(n);
            static_assert(std::is_same_v<typename tap::dsp::basic_real_fft<Sample>::engine,
                                         tap::dsp::detail::split_radix_rdft<Sample>>,
                          "basic_real_fft's declared engine is the split-radix port on this profile");
            for (const std::vector<Sample>& x : materials<Sample>(n)) {
                std::vector<Sample> via_class  = x;
                std::vector<Sample> via_engine = x;
                fft.forward_inplace(via_class.data());
                engine.forward_inplace(via_engine.data());
                ASSERT_EQ(std::memcmp(via_class.data(), via_engine.data(), n * sizeof(Sample)), 0)
                    << "forward N=" << n << ": basic_real_fft is not routed to split_radix_rdft byte for byte";

                std::vector<Sample> back_class  = via_engine;
                std::vector<Sample> back_engine = via_engine;
                fft.inverse_inplace(back_class.data());
                engine.inverse_inplace(back_engine.data());
                ASSERT_EQ(std::memcmp(back_class.data(), back_engine.data(), n * sizeof(Sample)), 0)
                    << "inverse N=" << n << ": basic_real_fft is not routed to split_radix_rdft byte for byte";
            }
        }
    }

    TEST(fft_routing, DoubleIsTheSplitRadixEngineByteForByte) {
        expect_routed_to_the_engine<double>();
    }

#if !defined(TAP_DSP_FFT_FLOAT_BACKEND)
    TEST(fft_routing, FloatIsTheSplitRadixEngineByteForByte) {
        expect_routed_to_the_engine<float>();
    }
#else
    // A backend build: basic_real_fft<float> is CMSIS or vDSP by design.
    // Pinned to the engine at float epsilon by tests/test_fft_backend.cpp.
    TEST(fft_routing, FloatIsABackendOnThisBuild) {
        static_assert(
            !std::is_same_v<tap::dsp::basic_real_fft<float>::engine, tap::dsp::detail::split_radix_rdft<float>>,
            "a backend define is active, so the float profile must not declare the port as its engine");
        SUCCEED();
    }
#endif

} // namespace
