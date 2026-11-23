// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Pins butterworth.h over iir_kernels.h's biquad cascade: the design's
// magnitude at DC, the cutoff and Nyquist, and its sections as the cookbook
// biquads at the Butterworth Qs; the stage against a double-precision direct
// recursion; chunking and channel independence bit for bit; the group delay,
// its nominal and its passband deviation (the IIR plan's I5); the tail after
// window_frames() bounded under the host's FTZ (I6), the window's values and
// its range at the cutoff bounds, and the limit cycle FTZ leaves some
// designs in, pinned as the honest limit; float's bound, where rounding
// stops threatening the poles; set_cutoff() as the constructor's design bit
// for bit, leaving the state alone, and the stage allocating nothing; the
// stage in a chain (in sequence, the nominal latency with a decimator, and
// flush equal to zero padding up to the recursive residue); that the stage
// never writes the floating-point environment (the host owns it); and the
// rejections. Every recursive test runs as a host would, FTZ set
// (iir_support.h's fixture); the same bounds without it are
// test_butterworth_subnormals.cpp's, in the tap_dsp_iir_subnormals
// executable (TAP_DSP_ALLOW_SUBNORMALS).
//
// This translation unit carries the heap meter's replacement of the global
// operator new and delete (tap/dsp/analysis/heap_meter.h) for its executable,
// tap_dsp_chain_stages (tests/CMakeLists.txt): exactly one per executable, and
// every allocation of the program, gtest's included, runs through it.

#define TAP_DSP_HEAP_METER_REPLACE_GLOBAL_NEW
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "support/iir_support.h"
#include "support/signals.h"
#include "tap/dsp/analysis/heap_meter.h"
#include "tap/dsp/butterworth.h"
#include "tap/dsp/chain.h"
#include "tap/dsp/decimate.h"
#include "tap/dsp/subnormals.h"

namespace {

    using tap::dsp::basic_butterworth_lowpass;
    using tap::dsp::iir_section;
    using tap::dsp::test::flush_to_zero_fixture;

    class butterworth : public flush_to_zero_fixture {};
    template <typename S>
    class butterworth_test : public flush_to_zero_fixture {};
    using floating_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(butterworth_test, floating_types, );

    template <std::size_t K>
    std::complex<double> response(const std::array<iir_section, K>& sec, double f) {
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

    TEST_F(butterworth, MagnitudeIsUnityAtDcHalfPowerAtCutoffZeroAtNyquist) {
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

    TEST_F(butterworth, SectionsRunInIncreasingPoleRadius) {
        const auto s = tap::dsp::design_butterworth_lowpass<8>(0.1);
        for (std::size_t i = 1; i < s.size(); ++i) {
            EXPECT_LT(s[i - 1].a2, s[i].a2); // |p|^2 = a2
        }
    }

    /// The cookbook low-pass biquad (Bristow-Johnson's "Audio EQ Cookbook")
    /// at f0 = fc cycles per sample and quality q, normalized by a0.
    iir_section cookbook_lowpass(double fc, double q) {
        const double w0    = 2.0 * std::numbers::pi * fc;
        const double cw    = std::cos(w0);
        const double alpha = std::sin(w0) / (2.0 * q);
        const double a0    = 1.0 + alpha;
        return {(1.0 - cw) / 2.0 / a0, (1.0 - cw) / a0, (1.0 - cw) / 2.0 / a0, -2.0 * cw / a0, (1.0 - alpha) / a0};
    }

    template <std::size_t N>
    double largest_coefficient_difference(const std::array<iir_section, N / 2>& ours, double fc,
                                          const std::array<double, N / 2>& q) {
        double worst = 0.0;
        for (std::size_t i = 0; i < N / 2; ++i) {
            const iir_section c = cookbook_lowpass(fc, q[i]);
            for (const double d :
                 {ours[i].b0 - c.b0, ours[i].b1 - c.b1, ours[i].b2 - c.b2, ours[i].a1 - c.a1, ours[i].a2 - c.a2}) {
                worst = std::max(worst, std::fabs(d));
            }
        }
        return worst;
    }

    // In exact arithmetic each section is the cookbook low-pass biquad at
    // Q_k = 1 / alpha_k (the file header), so a consumer's cookbook cascade at
    // the Butterworth Qs is this design: at the exact Qs within double
    // rounding, and at the Qs to eight decimals, as filter tables print them
    // (0.54119610 and 1.30656296 at order 4), within what those decimals cost.
    TEST_F(butterworth, SectionsAreTheCookbookBiquadsAtTheButterworthQs) {
        constexpr std::array<double, 2> q4_printed{0.54119610, 1.30656296};
        constexpr std::array<double, 4> q8_printed{0.50979558, 0.60134489, 0.89997622, 2.56291545};
        std::array<double, 2>           q4{};
        std::array<double, 4>           q8{};
        for (std::size_t i = 0; i < 2; ++i) { // processing order: lowest Q first
            q4[i] = 1.0 / (2.0 * std::sin(static_cast<double>(2 * (1 - i) + 1) * std::numbers::pi / 8.0));
            EXPECT_NEAR(q4[i], q4_printed[i], 5e-9);
        }
        for (std::size_t i = 0; i < 4; ++i) {
            q8[i] = 1.0 / (2.0 * std::sin(static_cast<double>(2 * (3 - i) + 1) * std::numbers::pi / 16.0));
            EXPECT_NEAR(q8[i], q8_printed[i], 5e-9);
        }
        double printed = 0.0;
        for (const double fc : {0.45 / 2, 0.45 / 4, 0.45 / 8, 0.45 / 16, 0.01, 0.25, 0.4}) {
            const auto s4 = tap::dsp::design_butterworth_lowpass<4>(fc);
            const auto s8 = tap::dsp::design_butterworth_lowpass<8>(fc);
            EXPECT_LT(largest_coefficient_difference<4>(s4, fc, q4), 1e-13) << fc;
            EXPECT_LT(largest_coefficient_difference<8>(s8, fc, q8), 1e-13) << fc;
            printed = std::max({printed, largest_coefficient_difference<4>(s4, fc, q4_printed),
                                largest_coefficient_difference<8>(s8, fc, q8_printed)});
        }
        std::printf("[ measured ] cookbook biquads at the eight-decimal Qs: coefficients within %.2g\n", printed);
        EXPECT_LT(printed, 1e-8);
    }

    TYPED_TEST(butterworth_test, TracksTheDoubleRecursion) {
        using sample_type = TypeParam;
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
        using sample_type                           = TypeParam;
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
        using sample_type                                = TypeParam;
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

    TEST_F(butterworth, GroupDelayIsThePhaseDerivative) {
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

    TEST_F(butterworth, PassbandDeviationIsTheLargestStrayUpToTheCutoff) {
        const basic_butterworth_lowpass<double, 4> f(0.1);
        double                                     worst = 0.0;
        for (int i = 0; i <= 4000; ++i) {
            const double fr = 0.1 * i / 4000.0;
            worst = std::max(worst, std::fabs(f.group_delay_output_frames(fr) - f.nominal_latency_output_frames()));
        }
        // The stage's 1,024-point grid finds the same peak to 1e-4 frames.
        EXPECT_NEAR(f.latency_deviation_output_frames(), worst, 1e-4);
        EXPECT_GT(f.latency_deviation_output_frames(), 0.0);
        // The header's numbers: 2.51 frames on a nominal of 4.02 at 0.1.
        EXPECT_NEAR(f.nominal_latency_output_frames(), 4.02, 0.005);
        EXPECT_NEAR(f.latency_deviation_output_frames(), 2.51, 0.005);
    }

    // The deviation is computed when it is read, from the sections and the
    // cutoff: the member is the free function of the stage's own design, the
    // same number on every read, and after set_cutoff() the new design's.
    TEST_F(butterworth, TheDeviationIsTheFreeFunctionOfTheDesignOnEveryRead) {
        basic_butterworth_lowpass<float, 8, 2> f(0.1);
        const double                           first = f.latency_deviation_output_frames();
        EXPECT_EQ(first, tap::dsp::butterworth_latency_deviation(f.sections(), f.cutoff()));
        EXPECT_EQ(first, f.latency_deviation_output_frames());
        f.set_cutoff(0.45 / 4.0);
        const basic_butterworth_lowpass<float, 8, 2> fresh(0.45 / 4.0);
        EXPECT_EQ(f.latency_deviation_output_frames(), fresh.latency_deviation_output_frames());
        EXPECT_EQ(f.latency_deviation_output_frames(),
                  tap::dsp::butterworth_latency_deviation(fresh.sections(), fresh.cutoff()));
        EXPECT_NE(f.latency_deviation_output_frames(), first);
    }

    // After window_frames() zeros the tail is the host's arithmetic's
    // (iir_kernels.h's honest limit): exact silence, or under FTZ a limit
    // cycle just above the smallest normal value for some designs. Whatever
    // it is, it is below 256 times the smallest normal value, from a
    // full-scale square drive and from one near the format's largest value.
    template <typename S, std::size_t N>
    void tail_is_bounded(double fc) {
        basic_butterworth_lowpass<S, N> f(fc);
        const S                         bound = S(256) * std::numeric_limits<S>::min();
        for (const S amp : {S(1), static_cast<S>(std::numeric_limits<S>::max() / S(1e6))}) {
            f.reset();
            auto x = square<S>(4000, 5, amp);
            f.process(x.data(), x.size(), x.data());
            std::vector<S> z(2 * f.window_frames(), S(0));
            f.process(z.data(), z.size(), z.data());
            for (std::size_t i = f.window_frames(); i < z.size(); ++i) {
                ASSERT_LE(std::fabs(z[i]), bound) << "fc " << fc << " frame " << i;
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
            tail_is_bounded<sample_type, 4>(fc);
            if (is_float) {
                tail_is_bounded<sample_type, 8>(fc);
            }
        }
    }

    TEST_F(butterworth, WindowIsADesignTimeConstant) {
        // Values recorded with this design (x86-64; a change is a change to
        // the design or the window formula).
        EXPECT_EQ((basic_butterworth_lowpass<float, 4>(0.1).window_frames()), 780u);
        EXPECT_EQ((basic_butterworth_lowpass<double, 4>(0.1).window_frames()), 6207u);
        EXPECT_EQ((basic_butterworth_lowpass<float, 8>(0.25).window_frames()), 906u);
    }

    // The cutoff bounds keep the window in a 32-bit size_t (the file header):
    // the largest, order 16 in double, is 232,346,074 frames at 1e-5 and
    // 232,710,082 at 1/2 - 1e-5 on x86-64 (a range here: libm's last bits
    // can move the window by a few frames on another target).
    TEST_F(butterworth, WindowStaysWithinA32BitSizeAtTheBounds) {
        constexpr double lo = tap::dsp::k_butterworth_min_cutoff<double>;
        for (const double fc : {lo, 0.5 - lo}) {
            const std::size_t w = basic_butterworth_lowpass<double, 16>(fc).window_frames();
            EXPECT_GT(w, std::size_t{230000000}) << fc;
            EXPECT_LT(w, std::size_t{235000000}) << fc;
        }
        constexpr double  lo32 = tap::dsp::k_butterworth_min_cutoff<float>;
        const std::size_t w32  = basic_butterworth_lowpass<float, 16>(lo32).window_frames();
        EXPECT_GT(w32, std::size_t{3000000});
        EXPECT_LT(w32, std::size_t{3100000});
    }

    // Float's bound (the file header): rounding a1 and a2 to float costs each
    // section's stability margin, 1 + a2 - |a1| (4 K^2 / a0 at the low end), at
    // most 9e-8. From the bound up every order keeps at least half of the
    // designed margin (0.906 of it at worst, measured); below it the rounding
    // can take all of it, and at 3e-5 it does for every order, a pole on the
    // unit circle (the honest limit the bound keeps out).
    template <std::size_t N>
    double float_margin_fraction(double fc) {
        double worst = std::numeric_limits<double>::infinity();
        for (const iir_section& s : tap::dsp::design_butterworth_lowpass<N>(fc)) {
            const double a1f = static_cast<double>(static_cast<float>(s.a1));
            const double a2f = static_cast<double>(static_cast<float>(s.a2));
            worst            = std::min(worst, (1.0 + a2f - std::fabs(a1f)) / (1.0 + s.a2 - std::fabs(s.a1)));
        }
        return worst;
    }

    TEST_F(butterworth, FloatRoundingKeepsThePolesInsideDownToTheFloatBound) {
        constexpr double lo = tap::dsp::k_butterworth_min_cutoff<float>;
        for (const double fc : {lo, 0.5 - lo, 2.0 * lo}) {
            EXPECT_GT(float_margin_fraction<2>(fc), 0.5) << fc;
            EXPECT_GT(float_margin_fraction<4>(fc), 0.5) << fc;
            EXPECT_GT(float_margin_fraction<8>(fc), 0.5) << fc;
            EXPECT_GT(float_margin_fraction<16>(fc), 0.5) << fc;
        }
        for (const double fc : {3e-5, 0.5 - 3e-5}) {
            EXPECT_LE(float_margin_fraction<2>(fc), 0.0) << fc;
            EXPECT_LE(float_margin_fraction<4>(fc), 0.0) << fc;
            EXPECT_LE(float_margin_fraction<8>(fc), 0.0) << fc;
            EXPECT_LE(float_margin_fraction<16>(fc), 0.0) << fc;
        }
    }

    // The honest limit, pinned on every leg (none skips). Where a flush bit
    // governs float, the host's FTZ leaves an order-8 design at fc = 0.45 in
    // a limit cycle above the smallest normal value (12.8 times it, measured
    // on x86-64) that never ends; where none does (the soft-float M4) the
    // residue stays below the smallest normal value, as it does without FTZ
    // anywhere (test_butterworth_subnormals.cpp).
    TEST_F(butterworth, UnderFlushToZeroSomeDesignsNeverFallSilent) {
        constexpr bool                      bit = tap::dsp::k_subnormal_control<float>;
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

    // ---- set_cutoff(): the constructor's design, in place ----

    template <typename S, std::size_t N>
    void expect_same_design(const basic_butterworth_lowpass<S, N>& a, const basic_butterworth_lowpass<S, N>& b) {
        EXPECT_EQ(a.cutoff(), b.cutoff());
        for (std::size_t i = 0; i < N / 2; ++i) {
            EXPECT_EQ(a.sections()[i].b0, b.sections()[i].b0) << "section " << i;
            EXPECT_EQ(a.sections()[i].b1, b.sections()[i].b1) << "section " << i;
            EXPECT_EQ(a.sections()[i].b2, b.sections()[i].b2) << "section " << i;
            EXPECT_EQ(a.sections()[i].a1, b.sections()[i].a1) << "section " << i;
            EXPECT_EQ(a.sections()[i].a2, b.sections()[i].a2) << "section " << i;
        }
        EXPECT_EQ(a.window_frames(), b.window_frames());
        EXPECT_EQ(a.nominal_latency_output_frames(), b.nominal_latency_output_frames());
        EXPECT_EQ(a.latency_deviation_output_frames(), b.latency_deviation_output_frames());
        EXPECT_EQ(a.group_delay_output_frames(0.3 * a.cutoff()), b.group_delay_output_frames(0.3 * b.cutoff()));
    }

    /// A stage at `from`, run on signal, then set_cutoff(to), against a
    /// fresh stage at `to`: the same design, and after reset() the same
    /// output bit for bit.
    template <typename S, std::size_t N>
    void set_cutoff_is_construction(double from, double to) {
        SCOPED_TRACE(testing::Message() << "order " << N << ", " << from << " -> " << to);
        const auto                      x = tap::dsp::test::random_signal<S>(1200, 17u, 0.8);
        basic_butterworth_lowpass<S, N> moved(from);
        std::vector<S>                  a(x.size()), b(x.size());
        moved.process(x.data(), x.size(), a.data());
        moved.set_cutoff(to);
        const basic_butterworth_lowpass<S, N> fresh_design(to);
        expect_same_design(moved, fresh_design);
        basic_butterworth_lowpass<S, N> fresh(to);
        moved.reset();
        moved.process(x.data(), x.size(), a.data());
        fresh.process(x.data(), x.size(), b.data());
        EXPECT_EQ(a, b);
    }

    TYPED_TEST(butterworth_test, SetCutoffIsTheConstructorsDesignBitForBit) {
        using sample_type   = TypeParam;
        constexpr double lo = tap::dsp::k_butterworth_min_cutoff<sample_type>;
        // An oversampling kernel's cutoffs (0.45 / os), both ways, and the bounds.
        for (const auto& [from, to] : std::array<std::array<double, 2>, 6>{
                 {{0.225, 0.1125}, {0.1125, 0.028125}, {0.028125, 0.225}, {0.05, 0.45}, {0.3, lo}, {lo, 0.5 - lo}}}) {
            set_cutoff_is_construction<sample_type, 4>(from, to);
            set_cutoff_is_construction<sample_type, 8>(from, to);
        }
    }

    // set_cutoff() leaves the state alone: to the same cutoff mid-stream it
    // changes nothing, bit for bit; to another, the recursion continues from
    // its state under the new coefficients (the double reference recursion
    // switching sections at the same frame, within TracksTheDoubleRecursion's
    // bounds), not from silence.
    TYPED_TEST(butterworth_test, SetCutoffLeavesTheStateAlone) {
        using sample_type       = TypeParam;
        constexpr std::size_t n = 2000, half = 1000;
        const auto            x = tap::dsp::test::random_signal<sample_type>(n, 23u, 0.9);

        basic_butterworth_lowpass<sample_type, 4> whole(0.1), split(0.1);
        std::vector<sample_type>                  a(n), b(n);
        whole.process(x.data(), n, a.data());
        split.process(x.data(), half, b.data());
        split.set_cutoff(0.1);
        split.process(x.data() + half, n - half, b.data() + half);
        EXPECT_EQ(a, b);

        basic_butterworth_lowpass<sample_type, 4> switched(0.1);
        std::vector<sample_type>                  y(n);
        switched.process(x.data(), half, y.data());
        switched.set_cutoff(0.03);
        EXPECT_FALSE(switched.is_silent());
        switched.process(x.data() + half, n - half, y.data() + half);
        const auto s0 = tap::dsp::design_butterworth_lowpass<4>(0.1);
        const auto s1 = tap::dsp::design_butterworth_lowpass<4>(0.03);
        double     st[2][2]{};
        double     worst = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const auto& sec = i < half ? s0 : s1;
            double      v   = static_cast<double>(x[i]);
            for (std::size_t k = 0; k < 2; ++k) {
                const auto&  c = sec[k];
                const double o = c.b0 * v + st[k][0];
                st[k][0]       = c.b1 * v - c.a1 * o + st[k][1];
                st[k][1]       = c.b2 * v - c.a2 * o;
                v              = o;
            }
            worst = std::max(worst, std::fabs(v - static_cast<double>(y[i])));
        }
        EXPECT_LT(worst, (std::is_same_v<sample_type, float> ? 2e-6 : 1e-14));
    }

    // The stage allocates nothing: construction, set_cutoff() (the design,
    // the window's matrix powers and the group-delay scan) at every order
    // the aliases name, in either format, up to the bounds, process() and
    // reset(), through the global operator new this executable replaces.
    template <typename S, std::size_t N>
    void run_without_allocating(std::vector<S>& buffer) {
        constexpr double                lo = tap::dsp::k_butterworth_min_cutoff<S>;
        basic_butterworth_lowpass<S, N> f(0.2);
        for (const double fc : {0.45 / 2, 0.45 / 4, 0.45 / 8, 0.45 / 16, lo, 0.5 - lo, 0.2}) {
            f.set_cutoff(fc);
            f.process(buffer.data(), buffer.size(), buffer.data());
        }
        f.reset();
    }

    TEST_F(butterworth, TheStageAllocatesNothing) {
        ASSERT_TRUE(tap::dsp::analysis::heap_meter::installed());
        auto       buffer   = tap::dsp::test::random_signal<float>(256, 3u, 0.5);
        auto       buffer64 = tap::dsp::test::random_signal<double>(256, 3u, 0.5);
        const auto s        = [&] {
            const tap::dsp::analysis::heap_meter meter;
            run_without_allocating<float, 4>(buffer);
            run_without_allocating<float, 8>(buffer);
            run_without_allocating<double, 4>(buffer64);
            run_without_allocating<double, 8>(buffer64);
            return meter.stats();
        }();
        EXPECT_EQ(s.allocations, 0u);
        EXPECT_EQ(s.total_bytes, 0u);
    }

#if !defined(NDEBUG) && defined(GTEST_HAS_DEATH_TEST) && GTEST_HAS_DEATH_TEST
    // The precondition itself: a debug assertion (TAP_EXPECTS), never a throw,
    // compiled in only where it can fire and be caught (a Debug configure on
    // a host with death tests, as fft_engine_debug's).
    TEST(butterworth_debug, SetCutoffOutsideTheRangeAssertsInADebugBuild) {
        basic_butterworth_lowpass<double, 4> f(0.1);
        basic_butterworth_lowpass<float, 4>  g(0.1);
        EXPECT_DEATH(f.set_cutoff(0.5), "");
        EXPECT_DEATH(f.set_cutoff(0.0), "");
        EXPECT_DEATH(f.set_cutoff(std::numeric_limits<double>::quiet_NaN()), "");
        EXPECT_DEATH(g.set_cutoff(5e-5), ""); // inside double's range, outside float's
        EXPECT_DEATH(static_cast<void>(tap::dsp::design_butterworth_lowpass<4>(0.6)), "");
        f.set_cutoff(tap::dsp::k_butterworth_min_cutoff<double>);
        g.set_cutoff(tap::dsp::k_butterworth_min_cutoff<float>);
        EXPECT_EQ(g.cutoff(), 1e-4) << "at the bound it does not fire";
    }
#endif

    // ---- In a chain (chain.h; SampleRateTap's IIR plan, I5 and I7) ----

    using tap::dsp::basic_decimator;
    using tap::dsp::chain;

    // Two stereo stages, orders 4 and 8, chain as the two run in sequence,
    // for any chunking.
    TYPED_TEST(butterworth_test, InAChainTheStagesRunInSequence) {
        using sample_type       = TypeParam;
        using bw4               = basic_butterworth_lowpass<sample_type, 4, 2>;
        using bw8               = basic_butterworth_lowpass<sample_type, 8, 2>;
        using two               = chain<bw4, bw8>;
        constexpr std::size_t n = 1000;
        const auto            x = tap::dsp::test::random_signal<sample_type>(2 * n, 21u, 0.7);
        static_assert(two::k_up == 1 && two::k_down == 1);
        two                      c(bw4(0.15), bw8(0.3));
        std::vector<sample_type> y(2 * n);
        std::size_t              o = 0, chunk = 1;
        while (o < n) {
            const std::size_t m = std::min(chunk, n - o);
            ASSERT_EQ(c.process(x.data() + 2 * o, m, y.data() + 2 * o), m);
            o += m;
            chunk = chunk * 3 % 101 + 1;
        }
        bw4                      first(0.15);
        bw8                      second(0.3);
        std::vector<sample_type> t(2 * n), r(2 * n);
        first.process(x.data(), n, t.data());
        second.process(t.data(), n, r.data());
        EXPECT_EQ(y, r);
        EXPECT_EQ(c.outputs_for(77), 77u);
        EXPECT_EQ(c.frames_needed(77), 77u);
    }

    TEST_F(butterworth, InAChainTheNominalLatencyComposesWithADecimator) {
        using bw  = basic_butterworth_lowpass<double, 4>;
        using dec = basic_decimator<double, 2>;
        chain<bw, dec> c(bw(0.2), dec{});
        const bw       alone(0.2);
        const double   dec_delay = static_cast<double>(c.stage<1>().latency_input_samples()) / 2.0;
        // The Butterworth runs at twice the chain's output rate: its frames count half.
        EXPECT_DOUBLE_EQ(c.nominal_latency_output_frames(), alone.nominal_latency_output_frames() / 2.0 + dec_delay);
        EXPECT_DOUBLE_EQ(c.latency_deviation_output_frames(), alone.latency_deviation_output_frames() / 2.0);
        for (const double f : {0.0, 0.05, 0.15, 0.3}) {
            EXPECT_DOUBLE_EQ(c.group_delay_output_frames(f), alone.group_delay_output_frames(f / 2.0) / 2.0 + dec_delay)
                << f;
        }
    }

    // Flush over a chain of the stage and a decimator: zero padding up to the
    // stage's residue, from several phases, under the host's FTZ (the
    // fixture); after it, that residue rather than exact silence (the file
    // header of butterworth.h). Without FTZ the residue is subnormal
    // (test_butterworth_subnormals.cpp), inside the same bound. An all-FIR
    // chain keeps the bit-for-bit equality (test_chain.cpp).
    TYPED_TEST(butterworth_test, InAChainFlushEqualsZeroPaddingUpToTheResidue) {
        using sample_type = TypeParam;
        using bw          = basic_butterworth_lowpass<sample_type, 4>;
        using dec         = basic_decimator<sample_type, 2>;
        using c2          = chain<bw, dec>;
        const auto x      = tap::dsp::test::random_signal<sample_type>(400, 8u, 0.9);
        for (std::size_t fed = 389; fed <= 400; fed += 3) {
            c2                       flushed(bw(0.12), dec{});
            std::vector<sample_type> scratch(flushed.outputs_for(fed));
            flushed.process(x.data(), fed, scratch.data());
            c2                       padded = flushed;
            const std::size_t        expect = flushed.flush_output_frames();
            std::vector<sample_type> a(expect + 4, sample_type(1));
            ASSERT_EQ(flushed.flush(a.data()), expect);
            const std::size_t        zeros = 2 * padded.window_frames();
            std::vector<sample_type> z(zeros, sample_type(0));
            std::vector<sample_type> b(padded.outputs_for(zeros));
            const std::size_t        made = padded.process(z.data(), zeros, b.data());
            ASSERT_GT(made, expect);
            // The stage's residue (below 256 times the smallest normal
            // value), through the decimator's gain (the sum of |h|, below 2):
            // flush and padding agree up to it (padding keeps feeding the
            // residue on, flush feeds the decimator pure zeros), and the
            // padded stream after the flush is at most it.
            const sample_type residue = sample_type(1024) * std::numeric_limits<sample_type>::min();
            for (std::size_t i = 0; i < expect; ++i) {
                ASSERT_LE(std::fabs(a[i] - b[i]), residue) << "fed " << fed << " output " << i;
            }
            for (std::size_t i = expect; i < made; ++i) {
                ASSERT_LE(std::fabs(b[i]), residue) << "fed " << fed << " output " << i;
            }
        }
    }

    template <typename T>
    concept has_window = requires(const T& t) { t.window_frames(); };
    static_assert(has_window<basic_butterworth_lowpass<float, 4>>);
    static_assert(tap::dsp::sync_stage<basic_butterworth_lowpass<float, 4, 2>>);
    static_assert(tap::dsp::sync_stage<basic_butterworth_lowpass<double, 8, 12>>);
    static_assert(tap::dsp::reports_nominal_latency<basic_butterworth_lowpass<double, 4>>);
    static_assert(std::is_same_v<tap::dsp::butterworth4<2>, basic_butterworth_lowpass<double, 4, 2>>);
    static_assert(std::is_same_v<tap::dsp::butterworth4_32<>, basic_butterworth_lowpass<float, 4, 1>>);
    static_assert(std::is_same_v<tap::dsp::butterworth8<>, basic_butterworth_lowpass<double, 8, 1>>);
    static_assert(std::is_same_v<tap::dsp::butterworth8_32<6>, basic_butterworth_lowpass<float, 8, 6>>);
    static_assert(noexcept(std::declval<basic_butterworth_lowpass<float, 8>&>().set_cutoff(0.1)));

    // The library leaves the floating-point environment to its host: it
    // reads the control register (a debug-build precondition) and never
    // writes it. Without FTZ: test_butterworth_subnormals.cpp.
    TEST_F(butterworth, NeverWritesTheFloatingPointEnvironment) {
        using tap::dsp::test::k_control_mask;
        const auto                          before = tap::dsp::floating_point_control() & k_control_mask;
        basic_butterworth_lowpass<float, 4> a(0.1);
        auto                                x = square<float>(512, 3);
        a.process(x.data(), x.size(), x.data());
        a.set_cutoff(0.2);
        a.reset();
        EXPECT_EQ(tap::dsp::floating_point_control() & k_control_mask, before);
    }

    TEST_F(butterworth, ValidateNamesTheCutoffRule) {
        using bw   = basic_butterworth_lowpass<double, 4>;
        using bw32 = basic_butterworth_lowpass<float, 4>;
        EXPECT_EQ(bw::validate(0.1), nullptr);
        EXPECT_EQ(bw::validate(1e-5), nullptr);
        EXPECT_EQ(bw::validate(0.5 - 1e-5), nullptr);
        EXPECT_NE(bw::validate(0.99e-5), nullptr);
        EXPECT_NE(bw::validate(0.0), nullptr);
        EXPECT_NE(bw::validate(0.5), nullptr);
        EXPECT_NE(bw::validate(-0.1), nullptr);
        EXPECT_NE(bw::validate(std::numeric_limits<double>::quiet_NaN()), nullptr);
        // Float's range is narrower (the file header: its rounding and the poles).
        EXPECT_EQ(bw32::validate(1e-4), nullptr);
        EXPECT_EQ(bw32::validate(0.5 - 1e-4), nullptr);
        EXPECT_NE(bw32::validate(5e-5), nullptr);
        EXPECT_NE(bw32::validate(0.49995), nullptr);
        static_assert(tap::dsp::butterworth_supports_cutoff<double>(5e-5));
        static_assert(!tap::dsp::butterworth_supports_cutoff<float>(5e-5));
#if TAP_DSP_HAS_EXCEPTIONS
        EXPECT_THROW(bw(0.5), std::invalid_argument);
        EXPECT_THROW(bw(0.0), std::invalid_argument);
        EXPECT_THROW(bw32(5e-5), std::invalid_argument);
#endif
    }

} // namespace
