// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Locks down tap::dsp::chain's contract (chain.h) on an existing primitive,
// as the rational engine's plan asks (M0): a chain of two decimators is the
// two run in sequence bit for bit, for any chunking, from a fresh state and
// after reset; outputs_for composes to exactly what process writes from
// every phase; frames_needed is the exact inverse; the latency is the exact
// rational the stages' delays compose to; flush drains every stage's tail and
// equals zero padding bit for bit, from every phase, in every sample format;
// and the concept admits basic_decimator and rejects a type missing a member.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/chain.h"
#include "tap/dsp/decimate.h"

namespace {

    using tap::dsp::basic_decimator;
    using tap::dsp::chain;
    using tap::dsp::decimate_profile;
    using tap::dsp::exact_ratio;
    using tap::dsp::sync_stage;

    static_assert(sync_stage<basic_decimator<float, 2>>);
    static_assert(sync_stage<basic_decimator<double, 3>>);
    static_assert(sync_stage<basic_decimator<std::int16_t, 6>>);

    struct not_a_stage {
        using sample                      = float;
        static constexpr std::size_t k_up = 1;
        // no k_down, no process
        void reset() noexcept {}
    };
    static_assert(!sync_stage<not_a_stage>);

    template <typename S>
    class chain_test : public ::testing::Test {};
    using sample_types = ::testing::Types<float, double, std::int16_t, std::int32_t>;
    TYPED_TEST_SUITE(chain_test, sample_types, );

    template <typename S>
    std::vector<S> test_signal(std::size_t n) {
        std::vector<S> x(n);
        for (std::size_t i = 0; i < n; ++i) {
            const double t = static_cast<double>(i);
            const double v = 0.4 * std::sin(2.0 * std::numbers::pi * 0.0113 * t)
                             + 0.2 * std::cos(2.0 * std::numbers::pi * 0.3 * t)
                             + 0.1 * std::sin(2.0 * std::numbers::pi * 0.47 * t);
            if constexpr (std::is_floating_point_v<S>) {
                x[i] = static_cast<S>(v);
            }
            else {
                x[i] = static_cast<S>(std::llround(v * (std::is_same_v<S, std::int16_t> ? 32767.0 : 2147483647.0)));
            }
        }
        return x;
    }

    // The reference: the two decimators run one after the other, whole.
    template <typename S>
    std::vector<S> reference_by_6(const std::vector<S>& x) {
        basic_decimator<S, 2> d2;
        basic_decimator<S, 3> d3;
        std::vector<S>        mid(d2.outputs_for(x.size()));
        const std::size_t     made = d2.process(x.data(), x.size(), mid.data());
        std::vector<S>        y(d3.outputs_for(made));
        const std::size_t     out = d3.process(mid.data(), made, y.data());
        y.resize(out);
        return y;
    }

    TYPED_TEST(chain_test, TwoDecimatorsInAChainAreTheTwoRunInSequence) {
        using sample                                                      = TypeParam;
        const auto                                                    x   = test_signal<sample>(5003);
        const auto                                                    ref = reference_by_6(x);
        chain<basic_decimator<sample, 2>, basic_decimator<sample, 3>> c(basic_decimator<sample, 2>{},
                                                                        basic_decimator<sample, 3>{});
        static_assert(decltype(c)::k_up == 1 && decltype(c)::k_down == 6);
        ASSERT_EQ(c.outputs_for(x.size()), ref.size());
        std::vector<sample> y(ref.size());
        ASSERT_EQ(c.process(x.data(), x.size(), y.data()), ref.size());
        EXPECT_TRUE(y == ref); // bit-exact, every sample type
    }

    TYPED_TEST(chain_test, ChunkingIsBitIdentical) {
        using sample   = TypeParam;
        const auto x   = test_signal<sample>(4001);
        const auto ref = reference_by_6(x);
        for (const std::size_t chunk : {1u, 5u, 63u, 64u, 65u, 700u, 4001u}) {
            chain<basic_decimator<sample, 2>, basic_decimator<sample, 3>> c(basic_decimator<sample, 2>{},
                                                                            basic_decimator<sample, 3>{});
            std::vector<sample>                                           y;
            y.reserve(ref.size());
            std::vector<sample> buf(ref.size() + 8);
            for (std::size_t pos = 0; pos < x.size(); pos += chunk) {
                const std::size_t n    = pos + chunk <= x.size() ? chunk : x.size() - pos;
                const std::size_t want = c.outputs_for(n);
                const std::size_t made = c.process(x.data() + pos, n, buf.data());
                ASSERT_EQ(made, want) << "chunk " << chunk << " at " << pos;
                y.insert(y.end(), buf.begin(), buf.begin() + static_cast<std::ptrdiff_t>(made));
            }
            EXPECT_TRUE(y == ref) << "chunk " << chunk;
        }
    }

    TEST(Chain, OutputsForIsExactFromEveryPhase) {
        chain<basic_decimator<float, 2>, basic_decimator<float, 3>> c(basic_decimator<float, 2>{},
                                                                      basic_decimator<float, 3>{});
        std::vector<float>                                          x(100, 0.25f), y(100);
        for (std::size_t phase = 0; phase < 12; ++phase) { // the chain's period is 6; walk two of them
            for (std::size_t n = 0; n <= 20; ++n) {
                chain<basic_decimator<float, 2>, basic_decimator<float, 3>> probe = c;
                EXPECT_EQ(probe.process(x.data(), n, y.data()), c.outputs_for(n)) << "phase " << phase << " n " << n;
            }
            c.process(x.data(), 1, y.data());
        }
    }

    TEST(Chain, FramesNeededIsTheExactInverseOfOutputsFor) {
        chain<basic_decimator<float, 2>, basic_decimator<float, 3>> c(basic_decimator<float, 2>{},
                                                                      basic_decimator<float, 3>{});
        std::vector<float>                                          x(64, 0.0f), y(64);
        for (std::size_t phase = 0; phase < 12; ++phase) {
            EXPECT_EQ(c.frames_needed(0), 0u);
            for (std::size_t k = 1; k <= 25; ++k) {
                const std::size_t n = c.frames_needed(k);
                EXPECT_GE(c.outputs_for(n), k) << "phase " << phase << " k " << k;
                EXPECT_LT(c.outputs_for(n - 1), k) << "phase " << phase << " k " << k;
            }
            c.process(x.data(), 1, y.data());
        }
    }

    TEST(Chain, LatencyIsTheExactRationalOfTheStagesDelays) {
        // by-2 economy: 81 taps, delay 40 input samples = 40/2 of its outputs;
        // by-3 economy: 121 taps, delay 60 of its inputs = 20 of its outputs.
        // At the chain's output rate: 40/2 * 1/3 + 60/3 = 20/3 + 20 = 80/3.
        chain<basic_decimator<float, 2>, basic_decimator<float, 3>> c(basic_decimator<float, 2>{},
                                                                      basic_decimator<float, 3>{});
        EXPECT_EQ(c.latency_output_frames(), (exact_ratio{80, 3}));
        EXPECT_NEAR(c.latency_seconds(16000.0), 80.0 / 3.0 / 16000.0, 1e-15);
        // Verified by impulse: the peak of the chain's response to a unit
        // impulse lands at output floor(80/3) or its neighbour, with the two
        // candidate samples straddling the exact delay.
        std::vector<float> x(400, 0.0f);
        x[0] = 1.0f;
        std::vector<float> y(c.outputs_for(x.size()));
        c.process(x.data(), x.size(), y.data());
        std::size_t peak = 0;
        for (std::size_t i = 1; i < y.size(); ++i) {
            if (y[i] > y[peak]) {
                peak = i;
            }
        }
        EXPECT_TRUE(peak == 26 || peak == 27) << peak;
    }

    TEST(Chain, SingleStageChainIsTheStage) {
        const auto                x = test_signal<float>(1000);
        basic_decimator<float, 3> d;
        std::vector<float>        ref(d.outputs_for(x.size()));
        d.process(x.data(), x.size(), ref.data());
        chain<basic_decimator<float, 3>> c(basic_decimator<float, 3>{});
        static_assert(decltype(c)::k_stages == 1);
        std::vector<float> y(c.outputs_for(x.size()));
        ASSERT_EQ(y.size(), ref.size());
        c.process(x.data(), x.size(), y.data());
        EXPECT_TRUE(y == ref);
        EXPECT_EQ(c.latency_output_frames(), (exact_ratio{20, 1}));
    }

    TEST(Chain, ResetReproducesBitExactly) {
        const auto                                                  x = test_signal<float>(777);
        chain<basic_decimator<float, 2>, basic_decimator<float, 3>> c(basic_decimator<float, 2>{},
                                                                      basic_decimator<float, 3>{});
        std::vector<float>                                          a(c.outputs_for(x.size()));
        c.process(x.data(), x.size(), a.data());
        std::vector<float> scratch(c.outputs_for(300));
        c.process(x.data(), 300, scratch.data()); // leave it mid-stream
        c.reset();
        std::vector<float> b(c.outputs_for(x.size()));
        ASSERT_EQ(a.size(), b.size());
        c.process(x.data(), x.size(), b.data());
        EXPECT_TRUE(a == b);
    }

    // flush() is zero padding: its output is a prefix of the padded stream's,
    // flush_output_frames() long, and the padded stream is silence beyond it.
    TYPED_TEST(chain_test, FlushEqualsZeroPaddingBitForBit) {
        using sample = TypeParam;
        using chain6 = chain<basic_decimator<sample, 2>, basic_decimator<sample, 3>>;
        const auto x = test_signal<sample>(601);
        for (std::size_t fed = 595; fed <= 601; ++fed) { // every phase of the by-6 chain
            chain6              flushed(basic_decimator<sample, 2>{}, basic_decimator<sample, 3>{});
            std::vector<sample> scratch(flushed.outputs_for(fed));
            flushed.process(x.data(), fed, scratch.data());
            chain6 padded = flushed;

            const std::size_t   expect = flushed.flush_output_frames();
            std::vector<sample> a(expect + 8, sample{1});
            ASSERT_EQ(flushed.flush(a.data()), expect);
            // Enough zeros to drain both stages twice over.
            const std::size_t zeros =
                2 * (padded.template stage<0>().window_frames() + 2 * padded.template stage<1>().window_frames());
            std::vector<sample> z(zeros, tap::dsp::sample_traits<sample>::silence());
            std::vector<sample> b(padded.outputs_for(zeros));
            const std::size_t   made = padded.process(z.data(), zeros, b.data());
            ASSERT_GT(made, expect);
            for (std::size_t i = 0; i < expect; ++i) {
                ASSERT_EQ(a[i], b[i]) << "fed " << fed << " output " << i;
            }
            for (std::size_t i = expect; i < made; ++i) {
                ASSERT_EQ(b[i], tap::dsp::sample_traits<sample>::silence()) << "fed " << fed << " output " << i;
            }
        }
    }

    TEST(Chain, FlushOutputFramesIsExactFromEveryPhaseAndTheTailIsDrained) {
        chain<basic_decimator<float, 2>, basic_decimator<float, 3>> c(basic_decimator<float, 2>{},
                                                                      basic_decimator<float, 3>{});
        // Each window is the stage's history length; the count is what
        // outputs_for composes over the two windows.
        const std::size_t w0 = c.stage<0>().window_frames();
        const std::size_t w1 = c.stage<1>().window_frames();
        EXPECT_EQ(w0, c.stage<0>().taps() - 1);
        EXPECT_EQ(w1, c.stage<1>().taps() - 1);
        EXPECT_EQ(c.flush_output_frames(), c.stage<1>().outputs_for(c.stage<0>().outputs_for(w0) + w1));
        const auto         x = test_signal<float>(64);
        std::vector<float> y(c.outputs_for(64) + c.flush_output_frames() + 1);
        for (std::size_t fed = 0; fed < 12; ++fed) {
            chain<basic_decimator<float, 2>, basic_decimator<float, 3>> probe = c;
            probe.process(x.data(), fed, y.data());
            const std::size_t expect = probe.flush_output_frames();
            EXPECT_EQ(probe.flush(y.data()), expect) << "fed " << fed;
            // Drained: the next window of zeros produces silence only.
            std::vector<float> z(w0 + 3 * w1, 0.0f);
            std::vector<float> tail(probe.outputs_for(z.size()));
            probe.process(z.data(), z.size(), tail.data());
            for (const float v : tail) {
                EXPECT_EQ(v, 0.0f);
            }
        }
    }

    TEST(Chain, ProfilesPassThroughToTheStages) {
        chain<basic_decimator<float, 2>, basic_decimator<float, 3>> c(
            basic_decimator<float, 2>{decimate_profile::transparent()}, basic_decimator<float, 3>{});
        EXPECT_EQ(c.stage<0>().taps(), 259u);
        EXPECT_EQ(c.stage<1>().taps(), 121u);
        EXPECT_EQ(c.channels(), 1u);
    }

    TEST(Chain, ExactRatioReducesAndComposes) {
        static_assert(exact_ratio{6, 4} == exact_ratio{3, 2});
        static_assert(exact_ratio{1, 2} * exact_ratio{1, 3} == exact_ratio{1, 6});
        static_assert(exact_ratio{1, 2} + exact_ratio{1, 3} == exact_ratio{5, 6});
        static_assert(exact_ratio{160, 147} * exact_ratio{1, 2} == exact_ratio{80, 147});
        EXPECT_DOUBLE_EQ((exact_ratio{80, 3}).value(), 80.0 / 3.0);
    }

} // namespace
