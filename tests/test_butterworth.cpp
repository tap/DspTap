// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Pins butterworth.h and the biquad kernel under it (biquad.h): the design's
// magnitude at DC, the cutoff and Nyquist; the stage against a double-
// precision direct recursion; chunking and channel independence bit for bit;
// the group delay, its nominal and its passband deviation (the IIR plan's
// I5); the tail after window_frames() bounded with the client's FTZ on and
// off (I6), and the limit cycle FTZ leaves some designs in, pinned as the
// honest limit; that the stage never touches the floating-point environment
// (the client owns it); and the rejections.

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <vector>

#include <gtest/gtest.h>

#include "support/flush_to_zero.h"
#include "support/signals.h"
#include "tap/dsp/butterworth.h"
#include "tap/dsp/chain.h"

namespace {

    using tap::dsp::basic_butterworth_lowpass;
    using tap::dsp::biquad_coefficients;
    using tap_test::scoped_flush_to_zero;

    template <typename S>
    class butterworth_test : public ::testing::Test {};
    using floating_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(butterworth_test, floating_types, );

    template <std::size_t K>
    std::complex<double> response(const std::array<biquad_coefficients, K>& sec, double f) {
        const auto           z1 = std::polar(1.0, -2.0 * std::numbers::pi * f);
        const auto           z2 = z1 * z1;
        std::complex<double> h  = 1.0;
        for (const auto& c : sec) {
            h *= (c.b0 + c.b1 * z1 + c.b2 * z2) / (1.0 + c.a1 * z1 + c.a2 * z2);
        }
        return h;
    }

    template <typename S>
    std::vector<S> square(std::size_t n, std::size_t half_period, S amp = S(1)) {
        std::vector<S> x(n);
        for (std::size_t i = 0; i < n; ++i) {
            x[i] = ((i / half_period) % 2) != 0 ? amp : -amp;
        }
        return x;
    }

    TEST(Butterworth, MagnitudeIsUnityAtDcHalfPowerAtCutoffZeroAtNyquist) {
        for (const double fc : {0.01, 0.05, 0.1, 0.25, 0.4, 0.49}) {
            const auto s4 = tap::dsp::design_butterworth_lowpass<4>(fc);
            EXPECT_NEAR(std::abs(response(s4, 0.0)), 1.0, 1e-12) << fc;
            EXPECT_NEAR(20.0 * std::log10(std::abs(response(s4, fc))), -3.0103, 1e-4) << fc;
            EXPECT_LT(std::abs(response(s4, 0.5)), 1e-12) << fc;
            const auto s8 = tap::dsp::design_butterworth_lowpass<8>(fc);
            EXPECT_NEAR(20.0 * std::log10(std::abs(response(s8, fc))), -3.0103, 1e-4) << fc;
            // Monotone: 2N log10 decay above the cutoff, a sample of it.
            EXPECT_LT(std::abs(response(s4, std::min(0.499, 1.5 * fc))), std::abs(response(s4, fc)));
        }
    }

    TEST(Butterworth, SectionsRunInIncreasingPoleRadius) {
        const auto s = tap::dsp::design_butterworth_lowpass<8>(0.1);
        for (std::size_t i = 1; i < s.size(); ++i) {
            EXPECT_LT(s[i - 1].a2, s[i].a2); // |p|^2 = a2
        }
    }

    TYPED_TEST(butterworth_test, TracksTheDoubleRecursion) {
        using sample_type = TypeParam;
        scoped_flush_to_zero                      ftz;
        basic_butterworth_lowpass<sample_type, 4> f(0.07);
        const auto                                x = tap::dsp::test::random_signal<double>(4000, 11u, 0.9);
        std::vector<sample_type>                  xs(x.size()), y(x.size());
        std::transform(x.begin(), x.end(), xs.begin(), [](double v) { return static_cast<sample_type>(v); });
        f.process(xs.data(), xs.size(), y.data());
        // The reference: the same sections, TDF II in double.
        const auto sec = tap::dsp::design_butterworth_lowpass<4>(0.07);
        double     st[2][2]{};
        double     worst = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            double v = static_cast<double>(xs[i]);
            for (std::size_t k = 0; k < 2; ++k) {
                const auto&  c = sec[k];
                const double o = c.b0 * v + st[k][0];
                st[k][0]       = c.b1 * v - c.a1 * o + st[k][1];
                st[k][1]       = c.b2 * v - c.a2 * o;
                v              = o;
            }
            worst = std::max(worst, std::fabs(v - static_cast<double>(y[i])));
        }
        // float: states and products in float (I14); double: bit-close.
        EXPECT_LT(worst, (std::is_same_v<sample_type, float> ? 2e-6 : 1e-14));
    }

    TYPED_TEST(butterworth_test, ChunkingIsBitIdentical) {
        using sample_type = TypeParam;
        scoped_flush_to_zero                      ftz;
        const auto                                x = tap::dsp::test::random_signal<sample_type>(3000, 5u, 0.8);
        basic_butterworth_lowpass<sample_type, 4> whole(0.11), pieces(0.11);
        std::vector<sample_type>                  a(x.size()), b(x.size());
        whole.process(x.data(), x.size(), a.data());
        std::size_t o = 0, chunk = 1;
        while (o < x.size()) {
            const std::size_t m = std::min(chunk, x.size() - o);
            pieces.process(x.data() + o, m, b.data() + o);
            o += m;
            chunk = chunk * 7 % 61 + 1;
        }
        EXPECT_EQ(a, b);
    }

    TYPED_TEST(butterworth_test, ChannelsAreIndependentAndInterleaved) {
        using sample_type = TypeParam;
        scoped_flush_to_zero                          ftz;
        constexpr std::size_t                         ch = 3, n = 1500;
        basic_butterworth_lowpass<sample_type, 4, ch> multi(0.2);
        std::vector<sample_type>                      in(ch * n), out(ch * n);
        std::array<std::vector<sample_type>, ch>      mono_in;
        for (std::size_t c = 0; c < ch; ++c) {
            mono_in[c] = tap::dsp::test::random_signal<sample_type>(n, static_cast<std::uint32_t>(c + 1), 0.7);
            for (std::size_t i = 0; i < n; ++i) {
                in[i * ch + c] = mono_in[c][i];
            }
        }
        multi.process(in.data(), n, out.data());
        for (std::size_t c = 0; c < ch; ++c) {
            basic_butterworth_lowpass<sample_type, 4, 1> mono(0.2);
            std::vector<sample_type>                     y(n);
            mono.process(mono_in[c].data(), n, y.data());
            for (std::size_t i = 0; i < n; ++i) {
                ASSERT_EQ(out[i * ch + c], y[i]) << "channel " << c << " frame " << i;
            }
        }
    }

    TEST(Butterworth, GroupDelayIsThePhaseDerivative) {
        const basic_butterworth_lowpass<double, 4> f(0.1);
        const auto                                 sec = f.sections();
        for (const double fr : {0.0, 0.02, 0.05, 0.09, 0.2}) {
            const double h  = 1e-6;
            const double f0 = std::max(fr - h, 0.0), f1 = fr + h;
            // -dphi/domega, with omega = 2 pi f, unwrapped locally.
            const double dphi = std::arg(response(sec, f1) / response(sec, f0));
            const double gd   = -dphi / (2.0 * std::numbers::pi * (f1 - f0));
            EXPECT_NEAR(f.group_delay_output_frames(fr), gd, 1e-4) << fr;
        }
        EXPECT_DOUBLE_EQ(f.nominal_latency_output_frames(), f.group_delay_output_frames(0.0));
    }

    TEST(Butterworth, PassbandDeviationIsTheLargestStrayUpToTheCutoff) {
        const basic_butterworth_lowpass<double, 4> f(0.1);
        double                                     worst = 0.0;
        for (int i = 0; i <= 4000; ++i) {
            const double fr = 0.1 * i / 4000.0;
            worst = std::max(worst, std::fabs(f.group_delay_output_frames(fr) - f.nominal_latency_output_frames()));
        }
        // The stage's 1,024-point grid finds the same peak to 1e-4 frames.
        EXPECT_NEAR(f.passband_group_delay_deviation_output_frames(), worst, 1e-4);
        EXPECT_GT(f.passband_group_delay_deviation_output_frames(), 0.0);
    }

    // After window_frames() zeros the tail is the client's arithmetic's
    // (biquad.h): exact silence, a subnormal residue without FTZ, or under
    // FTZ a limit cycle just above the smallest normal value for some
    // designs. Whatever it is, it is below 256 times the smallest normal
    // value, from a full-scale square drive and from one near the format's
    // largest value, with the client's FTZ on and off.
    template <typename S, std::size_t N>
    void tail_is_bounded(double fc, bool ftz_on) {
        scoped_flush_to_zero            mode(ftz_on);
        basic_butterworth_lowpass<S, N> f(fc);
        const S                         bound = S(256) * std::numeric_limits<S>::min();
        for (const S amp : {S(1), static_cast<S>(std::numeric_limits<S>::max() / S(1e6))}) {
            f.reset();
            auto x = square<S>(4000, 5, amp);
            f.process(x.data(), x.size(), x.data());
            std::vector<S> z(2 * f.window_frames(), S(0));
            f.process(z.data(), z.size(), z.data());
            for (std::size_t i = f.window_frames(); i < z.size(); ++i) {
                ASSERT_LE(std::fabs(z[i]), bound) << "fc " << fc << " ftz " << ftz_on << " frame " << i;
            }
        }
    }

    TYPED_TEST(butterworth_test, TheTailAfterTheWindowIsBounded) {
        using sample_type = TypeParam;
        // Double's windows are about 8 times float's (its range is), and
        // double is software on the M4, M4F and M33 legs: two cutoffs there.
        constexpr bool is_float = std::is_same_v<sample_type, float>;
        for (const double fc : {0.005, 0.02, 0.1, 0.25, 0.45}) {
            if (!is_float && (fc == 0.005 || fc == 0.02 || fc == 0.45)) {
                continue;
            }
            for (const bool ftz : {true, false}) {
                tail_is_bounded<sample_type, 4>(fc, ftz);
                if (is_float) {
                    tail_is_bounded<sample_type, 8>(fc, ftz);
                }
            }
        }
    }

    TEST(Butterworth, WindowIsADesignTimeConstant) {
        // Values recorded with this design (x86-64; a change is a change to
        // the design or the window formula).
        EXPECT_EQ((basic_butterworth_lowpass<float, 4>(0.1).window_frames()), 780u);
        EXPECT_EQ((basic_butterworth_lowpass<double, 4>(0.1).window_frames()), 6207u);
        EXPECT_EQ((basic_butterworth_lowpass<float, 8>(0.25).window_frames()), 906u);
    }

    // The honest limit, pinned both ways on every leg (none skips). Where a
    // flush bit governs float, the client's FTZ leaves an order-8 design at
    // fc = 0.45 in a limit cycle above the smallest normal value (12.8 times
    // it, measured on x86-64) that never ends; without FTZ (and on the
    // soft-float M4, which has no bit) the residue stays below the smallest
    // normal value.
    TEST(Butterworth, UnderFlushToZeroSomeDesignsNeverFallSilent) {
        constexpr bool                      bit = tap_test::has_flush_bit<float>();
        scoped_flush_to_zero                mode(bit);
        basic_butterworth_lowpass<float, 8> f(0.45);
        auto                                x = square<float>(4000, 5);
        f.process(x.data(), x.size(), x.data());
        std::vector<float> z(100000, 0.0f);
        f.process(z.data(), z.size(), z.data());
        float tail = 0.0f;
        for (std::size_t i = z.size() - 1000; i < z.size(); ++i) {
            tail = std::max(tail, std::fabs(z[i]));
        }
        if (bit) {
            EXPECT_FALSE(f.is_silent());
            EXPECT_GE(tail, std::numeric_limits<float>::min());
            EXPECT_LE(tail, 256.0f * std::numeric_limits<float>::min());
        }
        else {
            EXPECT_LT(tail, std::numeric_limits<float>::min());
        }
    }

    template <typename T>
    concept has_window = requires(const T& t) { t.window_frames(); };
    static_assert(has_window<basic_butterworth_lowpass<float, 4>>);
    static_assert(tap::dsp::sync_stage<basic_butterworth_lowpass<float, 4, 2>>);
    static_assert(tap::dsp::sync_stage<basic_butterworth_lowpass<double, 8, 12>>);

    // The library leaves the floating-point environment to its client.
    TEST(Butterworth, NeverTouchesTheFloatingPointEnvironment) {
        for (const bool on : {true, false}) {
            scoped_flush_to_zero                mode(on);
            const auto                          before = tap_test::fp_control_bits();
            basic_butterworth_lowpass<float, 4> a(0.1);
            auto                                x = square<float>(512, 3);
            a.process(x.data(), x.size(), x.data());
            a.reset();
            EXPECT_EQ(tap_test::fp_control_bits(), before);
        }
    }

    TEST(Butterworth, ValidateNamesTheCutoffRule) {
        using bw = basic_butterworth_lowpass<double, 4>;
        EXPECT_EQ(bw::validate(0.1), nullptr);
        EXPECT_NE(bw::validate(0.0), nullptr);
        EXPECT_NE(bw::validate(0.5), nullptr);
        EXPECT_NE(bw::validate(-0.1), nullptr);
        EXPECT_NE(bw::validate(std::numeric_limits<double>::quiet_NaN()), nullptr);
#if TAP_DSP_HAS_EXCEPTIONS
        EXPECT_THROW(bw(0.5), std::invalid_argument);
        EXPECT_THROW(bw(0.0), std::invalid_argument);
#endif
    }

} // namespace
