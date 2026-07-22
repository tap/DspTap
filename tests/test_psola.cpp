// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Locks down the tap::dsp::basic_psola contract: identity behavior at ratio 1
// (frequency and level preserved), shift accuracy across the ratio range with
// the period supplied by the caller, amplitude stability, clear() semantics,
// and float/double cross-precision agreement. The detector used as the pitch
// oracle is tap::dsp::yin, certified by its own battery.

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/psola.h"
#include "tap/dsp/yin.h"

namespace {

    constexpr double k_pi = 3.14159265358979323846;
    constexpr double k_sr = 48000.0;

    template <typename Sample>
    std::vector<Sample> run_sine(tap::dsp::basic_psola<Sample>& shifter, double freq, double ratio, double seconds) {
        const int           n      = static_cast<int>(seconds * k_sr);
        const Sample        period = static_cast<Sample>(k_sr / freq);
        std::vector<Sample> out(static_cast<size_t>(n));
        for (int t = 0; t < n; ++t) {
            const Sample x              = static_cast<Sample>(std::sin(2.0 * k_pi * freq * t / k_sr));
            out[static_cast<size_t>(t)] = shifter.process(x, period, static_cast<Sample>(ratio));
        }
        return out;
    }

    /// Band-limited sawtooth normalized to peak 1 — the harmonic-rich, voice-like
    /// material PSOLA is designed for (its spectral-envelope resampling needs
    /// harmonics to sample; see the header notes and the PureTone test below).
    template <typename Sample>
    std::vector<Sample> run_saw(tap::dsp::basic_psola<Sample>& shifter, double freq, double ratio, double seconds,
                                int harmonics = 20) {
        const int           n      = static_cast<int>(seconds * k_sr);
        const Sample        period = static_cast<Sample>(k_sr / freq);
        std::vector<double> wave(static_cast<size_t>(n), 0.0);
        double              peak = 0.0;
        for (int t = 0; t < n; ++t) {
            for (int h = 1; h <= harmonics; ++h) {
                wave[static_cast<size_t>(t)] += std::sin(2.0 * k_pi * freq * h * t / k_sr) / h;
            }
            peak = std::max(peak, std::abs(wave[static_cast<size_t>(t)]));
        }
        std::vector<Sample> out(static_cast<size_t>(n));
        for (int t = 0; t < n; ++t) {
            const Sample x              = static_cast<Sample>(wave[static_cast<size_t>(t)] / peak);
            out[static_cast<size_t>(t)] = shifter.process(x, period, static_cast<Sample>(ratio));
        }
        return out;
    }

    template <typename Sample>
    double measure_hz(const std::vector<Sample>& x) {
        const size_t        tau_min = static_cast<size_t>(k_sr / 2000.0);
        const size_t        tau_max = static_cast<size_t>(std::ceil(k_sr / 55.0));
        tap::dsp::yin       det(tau_max, tau_min, tau_max);
        std::vector<double> tail(det.frame_size());
        for (size_t i = 0; i < tail.size(); ++i) {
            tail[i] = static_cast<double>(x[x.size() - tail.size() + i]);
        }
        const auto r = det.analyze(tail.data());
        EXPECT_TRUE(r.voiced());
        return (r.period > 0.0) ? k_sr / r.period : 0.0;
    }

    template <typename Sample>
    double tail_peak(const std::vector<Sample>& x, size_t samples = 4800) {
        double peak = 0.0;
        for (size_t i = x.size() - samples; i < x.size(); ++i) {
            peak = std::max(peak, std::abs(static_cast<double>(x[i])));
        }
        return peak;
    }

    double cents(double f, double ref) {
        return 1200.0 * std::log2(f / ref);
    }

    template <typename Sample>
    class psola_test : public ::testing::Test {};

    using sample_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(psola_test, sample_types);

    TYPED_TEST(psola_test, GeometryAndLatency) {
        tap::dsp::basic_psola<TypeParam> shifter(900);
        EXPECT_EQ(shifter.max_period(), 900u);
        EXPECT_EQ(shifter.latency(), 2u * 900u + 2u);
    }

    TYPED_TEST(psola_test, IdentityRatioPreservesFrequencyAndLevel) {
        tap::dsp::basic_psola<TypeParam> shifter(900);
        const auto                       out = run_sine(shifter, 220.0, 1.0, 1.0);
        EXPECT_LT(std::abs(cents(measure_hz(out), 220.0)), 3.0);
        EXPECT_GT(tail_peak(out), 0.85);
        EXPECT_LT(tail_peak(out), 1.15);
    }

    TYPED_TEST(psola_test, ShiftAccuracyOnVoiceLikeMaterial) {
        for (const double ratio : {0.5, 0.8, 1.122462, 1.5, 2.0}) {
            tap::dsp::basic_psola<TypeParam> shifter(900);
            const auto                       out = run_saw(shifter, 150.0, ratio, 1.0);
            EXPECT_LT(std::abs(cents(measure_hz(out), 150.0 * ratio)), 8.0) << "ratio " << ratio;
            EXPECT_GT(tail_peak(out), 0.4) << "ratio " << ratio;
            // Downshifts concentrate the resampled envelope's harmonics and can
            // legitimately exceed the source peak; bound generously.
            EXPECT_LT(tail_peak(out), 2.2) << "ratio " << ratio;
        }
    }

    TYPED_TEST(psola_test, PureToneOctaveUpThinsOut) {
        // The documented envelope-sampling property: a pure sine shifted a full
        // octave leaves no harmonic at the envelope's only peak, so the output
        // all but vanishes. This is PSOLA being PSOLA (formant preservation),
        // not a defect — pinned so a change in this behavior is noticed.
        tap::dsp::basic_psola<TypeParam> shifter(900);
        const auto                       out = run_sine(shifter, 220.0, 2.0, 1.0);
        EXPECT_LT(tail_peak(out), 0.1);
    }

    TYPED_TEST(psola_test, OutputStaysFiniteOnNoiseLikePeriods) {
        tap::dsp::basic_psola<TypeParam> shifter(900);
        // Deliberately mismatched period (the caller's tracker can be wrong):
        // the shifter must stay bounded and finite regardless.
        for (int t = 0; t < 48000; ++t) {
            const TypeParam x = static_cast<TypeParam>(std::sin(2.0 * k_pi * 300.0 * t / k_sr));
            const TypeParam y = shifter.process(x, TypeParam(700), TypeParam(1.3));
            ASSERT_TRUE(std::isfinite(static_cast<double>(y)));
            ASSERT_LT(std::abs(static_cast<double>(y)), 4.0);
        }
    }

    TYPED_TEST(psola_test, ClearZerosTheState) {
        tap::dsp::basic_psola<TypeParam> shifter(900);
        run_sine(shifter, 220.0, 1.5, 0.25);
        shifter.clear();
        for (size_t i = 0; i < shifter.latency(); ++i) {
            EXPECT_EQ(shifter.process(TypeParam(0), TypeParam(218), TypeParam(1.5)), TypeParam(0));
        }
    }

    TEST(psola_cross_precision, FloatAgreesWithDoubleGoldenModel) {
        tap::dsp::psola   gold(900);
        tap::dsp::psola32 fast(900);
        const auto        od = run_saw(gold, 150.0, 1.5, 1.0);
        const auto        of = run_saw(fast, 150.0, 1.5, 1.0);
        const double      fd = measure_hz(od);
        const double      ff = measure_hz(of);
        EXPECT_NEAR(ff, fd, fd * 1e-3);
    }

} // namespace
