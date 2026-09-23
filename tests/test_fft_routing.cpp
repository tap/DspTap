// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// THE ROUTING PROOF (Stage 2b of docs/audit-fft-and-code-smells.md, Part 3;
// generalized at Stage 4): tap::dsp::basic_real_fft adds nothing between the
// caller and its engine — no scale, no copy that could round, no stale table
// — byte for byte, forward and inverse. Pinned for:
//
//   - basic_real_fft<double>, whose engine is detail::split_radix_rdft<double>
//     on every build;
//   - basic_real_fft<float, detail::split_radix_rdft<float>>, the split-radix
//     engine NAMED explicitly (Stage 4), on every build — including the M55
//     and macOS legs, where the float DEFAULT is CMSIS / vDSP and this
//     instantiation sits beside it in the same binary;
//   - basic_real_fft<float> (real_fft32), against whatever engine the build
//     selected as the default: where that is the split-radix engine the
//     check is the same memcmp; where it is an accelerated engine, the
//     class is pinned to that engine's own output byte for byte (the
//     engine's fidelity to the split-radix contract is test_fft_backend.cpp's
//     job, at float epsilon).
//
// Why this is a test and not a reading of fft.h. The bit-identity gate
// (tests/test_fft_parity_ooura.cpp) proves the ENGINE is the vendored C; it
// deliberately does not go through basic_real_fft. This file closes the
// remaining gap. Together the two files say what the Stage 2b PR claims:
// every consumer's output through basic_real_fft is the C's, bit for bit,
// on every profile that routes to the split-radix engine.
//
// memcmp, not EXPECT_EQ: == calls +0.0 and -0.0 equal and any NaN unequal to
// itself; the claim is about bits. Both sides receive identical input bits,
// and the inverse is judged on the engine's own forward output, so the two
// directions are pinned independently.

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
    // every code path the split-radix engine dispatches on: n = 4 (no real
    // post-pass, makect never called), 8 (the smallest post-pass), the
    // certified geometries 512 and 2048, and 4096 (the largest the emulated
    // legs run). An engine whose range excludes a size (CMSIS: 32 … 4096)
    // is checked at the sizes it supports.
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

    /// The class over `Engine` against the engine itself, byte for byte.
    template <typename Sample, typename Engine>
    void expect_routed_to_the_engine() {
        using fft = tap::dsp::basic_real_fft<Sample, Engine>;
        static_assert(std::is_same_v<typename fft::engine, Engine>, "the class declares the engine it was given");
        std::size_t checked = 0;
        for (const std::size_t n : k_sizes) {
            if (!fft::supports_size(n)) {
                continue;
            }
            fft    fft_under_test(n);
            Engine engine(n);
            for (const std::vector<Sample>& x : materials<Sample>(n)) {
                std::vector<Sample> via_class  = x;
                std::vector<Sample> via_engine = x;
                fft_under_test.forward_inplace(via_class.data());
                engine.forward_inplace(via_engine.data());
                ASSERT_EQ(std::memcmp(via_class.data(), via_engine.data(), n * sizeof(Sample)), 0)
                    << "forward N=" << n << ": basic_real_fft is not routed to its engine byte for byte";

                std::vector<Sample> back_class  = via_engine;
                std::vector<Sample> back_engine = via_engine;
                fft_under_test.inverse_inplace(back_class.data());
                engine.inverse_inplace(back_engine.data());
                ASSERT_EQ(std::memcmp(back_class.data(), back_engine.data(), n * sizeof(Sample)), 0)
                    << "inverse N=" << n << ": basic_real_fft is not routed to its engine byte for byte";
            }
            ++checked;
        }
        EXPECT_GE(checked, 3u) << "the engine's range must admit the certified geometries at least";
    }

    TEST(fft_routing, DoubleIsTheSplitRadixEngineByteForByte) {
        static_assert(std::is_same_v<tap::dsp::real_fft::engine, tap::dsp::detail::split_radix_rdft<double>>,
                      "double routes to the split-radix engine on every build");
        expect_routed_to_the_engine<double, tap::dsp::detail::split_radix_rdft<double>>();
    }

    // On every build, the M55 and macOS legs included: the split-radix
    // engine named explicitly is the engine, byte for byte.
    TEST(fft_routing, FloatNamedSplitRadixIsTheEngineByteForByte) {
        expect_routed_to_the_engine<float, tap::dsp::detail::split_radix_rdft<float>>();
    }

    // The float default: the split-radix engine where no accelerated engine
    // was selected (the Stage 2b claim, unchanged), otherwise the selected
    // engine's own output — and in either case the class adds nothing.
    TEST(fft_routing, FloatDefaultIsTheSelectedEngineByteForByte) {
        using selected = tap::dsp::default_real_fft_engine_t<float>;
        static_assert(std::is_same_v<tap::dsp::real_fft32::engine, selected>);
#if defined(TAP_DSP_FFT_CMSIS) || defined(TAP_DSP_FFT_ACCELERATE)
        static_assert(!std::is_same_v<selected, tap::dsp::detail::split_radix_rdft<float>>,
                      "a backend define is active, so the float default must not be the split-radix engine");
#else
        static_assert(std::is_same_v<selected, tap::dsp::detail::split_radix_rdft<float>>,
                      "no backend define is active, so the float default is the split-radix engine");
#endif
        expect_routed_to_the_engine<float, selected>();
    }

} // namespace
