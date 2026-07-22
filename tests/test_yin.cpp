// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Locks down the tap::dsp::basic_yin contract: sub-sample period accuracy on
// pure and harmonic-rich (octave-trap) material, unvoiced rejection of noise
// and silence, threshold semantics, and float/double cross-precision
// agreement. The TapTools tune kernel's correctness later depends on every one
// of these staying exactly as documented in yin.h.

#include <cmath>
#include <numbers>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/yin.h"

namespace {

    constexpr double k_sr      = 48000.0;
    constexpr size_t k_window  = 800;
    constexpr size_t k_tau_min = 20;  // 2400 Hz
    constexpr size_t k_tau_max = 800; // 60 Hz

    template <typename Sample>
    std::vector<Sample> sine(double freq, size_t n, double amp = 1.0, double phase = 0.0) {
        std::vector<Sample> x(n);
        for (size_t i = 0; i < n; ++i) {
            x[i] = static_cast<Sample>(
                amp * std::sin(phase + 2.0 * std::numbers::pi * freq * static_cast<double>(i) / k_sr));
        }
        return x;
    }

    // Band-limited sawtooth: the classic octave trap for naive autocorrelation,
    // which happily locks to 2x the period on harmonic-rich material.
    template <typename Sample>
    std::vector<Sample> saw(double freq, size_t n, int harmonics) {
        std::vector<Sample> x(n, Sample(0));
        for (int h = 1; h <= harmonics; ++h) {
            for (size_t i = 0; i < n; ++i) {
                const double w = 2.0 * std::numbers::pi * freq * static_cast<double>(h) * static_cast<double>(i) / k_sr;
                x[i] += static_cast<Sample>(std::sin(w) / static_cast<double>(h));
            }
        }
        return x;
    }

    template <typename Sample>
    std::vector<Sample> noise(size_t n, unsigned seed, double amp = 1.0) {
        std::mt19937                           gen(seed);
        std::uniform_real_distribution<double> dist(-amp, amp);
        std::vector<Sample>                    x(n);
        for (auto& v : x) {
            v = static_cast<Sample>(dist(gen));
        }
        return x;
    }

    double cents_error(double detected_period, double true_freq) {
        const double detected_freq = k_sr / detected_period;
        return 1200.0 * std::log2(detected_freq / true_freq);
    }

    // Sub-sample interpolation should land well under a cent on clean material
    // in double; the float profile accumulates the difference function in
    // float32, so it gets a wider (but still musically negligible) budget.
    template <typename Sample>
    constexpr double k_cents_tolerance = 0.0;
    template <>
    constexpr double k_cents_tolerance<double> = 1.0;
    template <>
    constexpr double k_cents_tolerance<float> = 3.0;

    template <typename Sample>
    class yin_test : public ::testing::Test {};

    using sample_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(yin_test, sample_types);

    TYPED_TEST(yin_test, GeometryAndDefaults) {
        tap::dsp::basic_yin<TypeParam> det(k_window, k_tau_min, k_tau_max);
        EXPECT_EQ(det.window(), k_window);
        EXPECT_EQ(det.tau_min(), k_tau_min);
        EXPECT_EQ(det.tau_max(), k_tau_max);
        EXPECT_EQ(det.frame_size(), k_window + k_tau_max);
        EXPECT_NEAR(det.threshold(), 0.1, 1e-7);
    }

    TYPED_TEST(yin_test, SinePeriodAccuracyAcrossRange) {
        tap::dsp::basic_yin<TypeParam> det(k_window, k_tau_min, k_tau_max);

        // Includes frequencies with deliberately non-integer periods (439.7,
        // 441.3) to exercise the parabolic sub-sample refinement.
        for (const double freq : {82.4, 110.0, 146.8, 220.0, 261.6, 439.7, 440.0, 441.3, 880.0, 987.8}) {
            const auto x = sine<TypeParam>(freq, det.frame_size());
            const auto r = det.analyze(x.data());
            ASSERT_TRUE(r.voiced()) << freq << " Hz";
            EXPECT_LT(std::abs(cents_error(static_cast<double>(r.period), freq)), k_cents_tolerance<TypeParam>)
                << freq << " Hz -> period " << r.period;
            EXPECT_LT(r.aperiodicity, TypeParam(0.02)) << freq << " Hz";
        }
    }

    TYPED_TEST(yin_test, SawtoothFindsFundamentalNotOctave) {
        tap::dsp::basic_yin<TypeParam> det(k_window, k_tau_min, k_tau_max);

        for (const double freq : {110.0, 220.0, 440.0}) {
            const auto x = saw<TypeParam>(freq, det.frame_size(), 20);
            const auto r = det.analyze(x.data());
            ASSERT_TRUE(r.voiced()) << freq << " Hz";
            // Within tolerance of the fundamental — and therefore nowhere near
            // the half/double-period octave errors (+-1200 cents).
            EXPECT_LT(std::abs(cents_error(static_cast<double>(r.period), freq)), k_cents_tolerance<TypeParam>)
                << freq << " Hz -> period " << r.period;
        }
    }

    TYPED_TEST(yin_test, NoiseAndSilenceAreUnvoiced) {
        tap::dsp::basic_yin<TypeParam> det(k_window, k_tau_min, k_tau_max);

        const auto n  = noise<TypeParam>(det.frame_size(), 42);
        const auto rn = det.analyze(n.data());
        EXPECT_FALSE(rn.voiced());
        EXPECT_GT(rn.aperiodicity, det.threshold());

        const std::vector<TypeParam> zeros(det.frame_size(), TypeParam(0));
        const auto                   rz = det.analyze(zeros.data());
        EXPECT_FALSE(rz.voiced());
    }

    TYPED_TEST(yin_test, ThresholdGatesVoicing) {
        tap::dsp::basic_yin<TypeParam> det(k_window, k_tau_min, k_tau_max);

        // A deliberately dirty sine: periodic enough to pass the default
        // threshold, dirty enough that a tightened threshold rejects it.
        auto       x = sine<TypeParam>(220.0, det.frame_size());
        const auto n = noise<TypeParam>(det.frame_size(), 7, 0.25);
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] += n[i];
        }

        const auto voiced = det.analyze(x.data());
        ASSERT_TRUE(voiced.voiced());
        EXPECT_GT(voiced.aperiodicity, TypeParam(0));

        det.set_threshold(voiced.aperiodicity * TypeParam(0.5));
        const auto rejected = det.analyze(x.data());
        EXPECT_FALSE(rejected.voiced());
    }

    TEST(yin_cross_precision, FloatAgreesWithDoubleGoldenModel) {
        tap::dsp::yin   gold(k_window, k_tau_min, k_tau_max);
        tap::dsp::yin32 fast(k_window, k_tau_min, k_tau_max);

        for (const double freq : {110.0, 261.6, 439.7, 880.0}) {
            const auto xd = sine<double>(freq, gold.frame_size());
            const auto xf = sine<float>(freq, fast.frame_size());
            const auto rd = gold.analyze(xd.data());
            const auto rf = fast.analyze(xf.data());
            ASSERT_TRUE(rd.voiced());
            ASSERT_TRUE(rf.voiced());
            EXPECT_NEAR(static_cast<double>(rf.period), rd.period, rd.period * 1e-3) << freq << " Hz";
        }
    }

} // namespace
