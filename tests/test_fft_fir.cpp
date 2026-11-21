// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Pins fft_fir.h: the output is the direct convolution delayed by exactly B
// frames (to the FFT's accuracy), for any chunking bit for bit, per channel
// independently through the interleaved frames; the exact latency for
// linear-phase taps; silence after window_frames() zeros from every phase of
// a block; and the construction rules.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <gtest/gtest.h>

#include "support/signals.h"
#include "tap/dsp/chain.h"
#include "tap/dsp/fft_fir.h"
#include "tap/dsp/kaiser.h"

namespace {

    using tap::dsp::basic_fft_fir;
    using tap::dsp::exact_ratio;

    template <typename S>
    class fft_fir_test : public ::testing::Test {};
    using floating_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(fft_fir_test, floating_types, );

    std::vector<double> lowpass(std::size_t taps, double cutoff = 0.3) {
        std::vector<double> h(taps);
        tap::dsp::design_prototype(h, 1, cutoff, tap::dsp::kaiser_beta(80.0));
        return h;
    }

    template <typename S, std::size_t Ch>
    void matches_direct(std::size_t taps, std::size_t block) {
        const auto            h = lowpass(taps);
        basic_fft_fir<S, Ch>  f(h, block);
        const std::size_t     b  = f.block_delay_frames();
        constexpr std::size_t n  = 1500;
        const auto            xd = tap::dsp::test::random_signal<double>(n * Ch, 9u, 0.9);
        std::vector<S>        x(xd.begin(), xd.end()), y(n * Ch);
        f.process(x.data(), n, y.data());
        double worst = 0.0;
        for (std::size_t t = 0; t < n; ++t) {
            for (std::size_t c = 0; c < Ch; ++c) {
                double d = 0.0;
                for (std::size_t k = 0; k < taps; ++k) {
                    if (t >= b + k) {
                        d += h[k] * static_cast<double>(x[(t - b - k) * Ch + c]);
                    }
                }
                worst = std::max(worst, std::fabs(d - static_cast<double>(y[t * Ch + c])));
            }
        }
        EXPECT_LT(worst, (std::is_same_v<S, float> ? 1e-6 : 1e-14)) << taps << " taps, block " << b;
        for (std::size_t t = 0; t < b * Ch; ++t) {
            ASSERT_EQ(y[t], S(0)); // the first block is silence
        }
    }

    TYPED_TEST(fft_fir_test, MatchesDirectConvolution) {
        using sample_type = TypeParam;
        matches_direct<sample_type, 1>(255, 0);
        matches_direct<sample_type, 1>(255, 1024);
        matches_direct<sample_type, 2>(64, 0);
        matches_direct<sample_type, 6>(33, 32);
        matches_direct<sample_type, 1>(1, 0);
        matches_direct<sample_type, 1>(17, 16);
    }

    TYPED_TEST(fft_fir_test, ChunkingIsBitIdentical) {
        using sample_type               = TypeParam;
        const auto                    h = lowpass(101);
        const auto                    x = tap::dsp::test::random_signal<sample_type>(2 * 2500, 3u, 0.8);
        basic_fft_fir<sample_type, 2> whole(h), pieces(h);
        std::vector<sample_type>      a(x.size()), b(x.size());
        whole.process(x.data(), 2500, a.data());
        std::size_t o = 0, chunk = 1;
        while (o < 2500) {
            const std::size_t m = std::min(chunk, std::size_t{2500} - o);
            pieces.process(x.data() + 2 * o, m, b.data() + 2 * o);
            o += m;
            chunk = chunk * 5 % 173 + 1;
        }
        EXPECT_EQ(a, b);
    }

    TYPED_TEST(fft_fir_test, InPlaceProcessingIsTheSame) {
        using sample_type            = TypeParam;
        const auto                 h = lowpass(64);
        auto                       x = tap::dsp::test::random_signal<sample_type>(700, 4u, 0.8);
        basic_fft_fir<sample_type> a(h), b(h);
        std::vector<sample_type>   y(x.size());
        a.process(x.data(), x.size(), y.data());
        b.process(x.data(), x.size(), x.data());
        EXPECT_EQ(x, y);
    }

    TYPED_TEST(fft_fir_test, ChannelsAreIndependentAndInterleaved) {
        using sample_type                 = TypeParam;
        constexpr std::size_t          ch = 3, n = 900;
        const auto                     h = lowpass(48);
        basic_fft_fir<sample_type, ch> multi(h);
        std::vector<sample_type>       in(ch * n), out(ch * n);
        for (std::size_t c = 0; c < ch; ++c) {
            const auto m = tap::dsp::test::random_signal<sample_type>(n, static_cast<std::uint32_t>(c + 7), 0.6);
            for (std::size_t i = 0; i < n; ++i) {
                in[i * ch + c] = m[i];
            }
        }
        multi.process(in.data(), n, out.data());
        for (std::size_t c = 0; c < ch; ++c) {
            basic_fft_fir<sample_type, 1> mono(h);
            std::vector<sample_type>      mi(n), mo(n);
            for (std::size_t i = 0; i < n; ++i) {
                mi[i] = in[i * ch + c];
            }
            mono.process(mi.data(), n, mo.data());
            for (std::size_t i = 0; i < n; ++i) {
                ASSERT_EQ(out[i * ch + c], mo[i]) << "channel " << c << " frame " << i;
            }
        }
    }

    TEST(FftFir, LatencyIsExactForLinearPhaseTaps) {
        const basic_fft_fir<double> odd(lowpass(255));
        EXPECT_TRUE(odd.linear_phase());
        EXPECT_EQ(odd.latency_output_frames(), (exact_ratio{256 + 127, 1}));
        const basic_fft_fir<double> even(lowpass(64));
        EXPECT_EQ(even.latency_output_frames(), (exact_ratio{2 * 64 + 63, 2}));
        std::vector<double> anti{1.0, 0.5, 0.0, -0.5, -1.0};
        EXPECT_TRUE(basic_fft_fir<double>(anti).linear_phase());
        std::vector<double>         skew{1.0, 0.5, 0.25};
        const basic_fft_fir<double> s(skew);
        EXPECT_FALSE(s.linear_phase());
        EXPECT_EQ(s.block_delay_frames(), 16u);
    }

    TYPED_TEST(fft_fir_test, WindowDrainsToSilenceFromEveryPhase) {
        using sample_type            = TypeParam;
        const auto                 h = lowpass(31);
        basic_fft_fir<sample_type> f(h, 32);
        for (std::size_t p = 0; p < 32; p += 5) {
            f.reset();
            auto                     x = tap::dsp::test::random_signal<sample_type>(100 + p, 2u, 1.0);
            std::vector<sample_type> y(x.size());
            f.process(x.data(), x.size(), y.data());
            std::vector<sample_type> z(f.window_frames(), sample_type(0));
            f.process(z.data(), z.size(), z.data());
            std::vector<sample_type> more(200, sample_type(0));
            f.process(more.data(), more.size(), more.data());
            for (const sample_type v : more) {
                ASSERT_EQ(v, sample_type(0)) << "phase " << p;
            }
        }
        EXPECT_EQ(f.window_frames(), 96u);
    }

    TEST(FftFir, ValidateNamesEachRule) {
        using fir = basic_fft_fir<double>;
        EXPECT_EQ(fir::validate(255, 256), nullptr);
        EXPECT_NE(fir::validate(0, 16), nullptr);
        EXPECT_NE(fir::validate(255, 255), nullptr); // not a power of two
        EXPECT_NE(fir::validate(255, 128), nullptr); // below taps - 1
        EXPECT_EQ(fir::default_block(255), 256u);
        EXPECT_EQ(fir::default_block(1), 16u);
        EXPECT_EQ(fir::default_block(257), 256u);
#if TAP_DSP_HAS_EXCEPTIONS
        const auto h = lowpass(255);
        EXPECT_THROW(fir(h, 128), std::invalid_argument);
#endif
    }

    static_assert(tap::dsp::sync_stage<basic_fft_fir<float, 6>>);
    static_assert(tap::dsp::sync_stage<basic_fft_fir<double, 12>>);

} // namespace
