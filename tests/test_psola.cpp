// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Locks down the tap::dsp::basic_psola contract: identity behavior at ratio 1
// (frequency and level preserved), shift accuracy across the ratio range with
// the period supplied by the caller, amplitude stability, clear() semantics,
// and float/double cross-precision agreement. The detector used as the pitch
// oracle is tap::dsp::yin, certified by its own battery.

#include <cmath>
#include <cstdint>
#include <numbers>
#include <vector>

#include <gtest/gtest.h>

#include "support/pitch.h"
#include "tap/dsp/psola.h"
#include "tap/dsp/yin.h"

namespace {

    using tap::dsp::test::cents;
    using tap::dsp::test::measure_hz;
    using tap::dsp::test::run_sine;
    using tap::dsp::test::tail_peak;

    constexpr double k_sr = 48000.0;

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
                wave[static_cast<size_t>(t)] += std::sin(2.0 * std::numbers::pi * freq * h * t / k_sr) / h;
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
        const auto                       out = run_sine(shifter, 220.0, 1.0, 1.0, k_sr);
        EXPECT_LT(std::abs(cents(measure_hz(out, k_sr), 220.0)), 3.0);
        EXPECT_GT(tail_peak(out), 0.85);
        EXPECT_LT(tail_peak(out), 1.15);
    }

    TYPED_TEST(psola_test, ShiftAccuracyOnVoiceLikeMaterial) {
        for (const double ratio : {0.5, 0.8, 1.122462, 1.5, 2.0}) {
            tap::dsp::basic_psola<TypeParam> shifter(900);
            const auto                       out = run_saw(shifter, 150.0, ratio, 1.0);
            EXPECT_LT(std::abs(cents(measure_hz(out, k_sr), 150.0 * ratio)), 8.0) << "ratio " << ratio;
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
        const auto                       out = run_sine(shifter, 220.0, 2.0, 1.0, k_sr);
        EXPECT_LT(tail_peak(out), 0.1);
    }

    TYPED_TEST(psola_test, OutputStaysFiniteOnNoiseLikePeriods) {
        tap::dsp::basic_psola<TypeParam> shifter(900);
        // Deliberately mismatched period (the caller's tracker can be wrong):
        // the shifter must stay bounded and finite regardless.
        for (int t = 0; t < 48000; ++t) {
            const TypeParam x = static_cast<TypeParam>(std::sin(2.0 * std::numbers::pi * 300.0 * t / k_sr));
            const TypeParam y = shifter.process(x, TypeParam(700), TypeParam(1.3));
            ASSERT_TRUE(std::isfinite(static_cast<double>(y)));
            ASSERT_LT(std::abs(static_cast<double>(y)), 4.0);
        }
    }

    TYPED_TEST(psola_test, ClearZerosTheState) {
        tap::dsp::basic_psola<TypeParam> shifter(900);
        run_sine(shifter, 220.0, 1.5, 0.25, k_sr);
        shifter.clear();
        for (size_t i = 0; i < shifter.latency(); ++i) {
            EXPECT_EQ(shifter.process(TypeParam(0), TypeParam(218), TypeParam(1.5)), TypeParam(0));
        }
    }

    TYPED_TEST(psola_test, ClockWrapIsSeamlessAndCountersDoNotOverflow) {
        // The run-time contract is a bound, not an observable: the 32-bit sample
        // clock stays below 2 * clock_wrap(), a multiple of the ring size, and
        // the wrap changes nothing but the magnitude of the fractional mark
        // positions. Reference: a shifter run from zero. Subject: the same
        // shifter with its clock advanced past 2^31 elapsed samples through the
        // documented O(1) seam (the count a 48 kHz Cortex-M or Windows build
        // reaches after 12.4 h), then run for more than one wrap period so the
        // wrap fires inside the run. Ratio 1.5 makes the synthesis step t / r
        // inexact, so the wrap's magnitude change is exercised; the outputs agree
        // to the rounding of those positions (measured 2e-9 against a 0.35 peak),
        // five orders below any dropped or doubled grain.
        using shifter_t = tap::dsp::basic_psola<TypeParam>;
        shifter_t ref(900);
        shifter_t sub(900);
        ASSERT_LE(sub.clock_wrap(), size_t{1} << 18);    // the documented span
        ASSERT_EQ(sub.clock_wrap() % (4 * 900 + 8), 0u); // a multiple of the ring size

        const TypeParam period = static_cast<TypeParam>(k_sr / 150.0);
        const int       warm   = 12000;            // seam applied past the warm-up, mid-stream
        const int       run    = (1 << 18) + 8192; // > clock_wrap(), so the wrap fires inside the run
        double          worst  = 0.0;
        for (int t = 0; t < warm + run; ++t) {
            if (t == warm) {
                sub.advance_clock_for_testing((std::uint64_t{1} << 31) + 12345);
            }
            double x = 0.0;
            for (int h = 1; h <= 20; ++h) {
                x += std::sin(2.0 * std::numbers::pi * 150.0 * h * t / k_sr) / h;
            }
            const TypeParam xs = static_cast<TypeParam>(x / 3.6);
            const double    a  = static_cast<double>(ref.process(xs, period, TypeParam(1.5)));
            const double    b  = static_cast<double>(sub.process(xs, period, TypeParam(1.5)));
            ASSERT_TRUE(std::isfinite(b)) << "t " << t;
            worst = std::max(worst, std::abs(a - b));
        }
        EXPECT_LT(worst, 1e-6);
    }

    TEST(psola_cross_precision, FloatAgreesWithDoubleGoldenModel) {
        tap::dsp::psola   gold(900);
        tap::dsp::psola32 fast(900);
        const auto        od = run_saw(gold, 150.0, 1.5, 1.0);
        const auto        of = run_saw(fast, 150.0, 1.5, 1.0);
        const double      fd = measure_hz(od, k_sr);
        const double      ff = measure_hz(of, k_sr);
        EXPECT_NEAR(ff, fd, fd * 1e-3);
    }

} // namespace
