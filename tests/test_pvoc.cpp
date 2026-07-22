// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Locks down the tap::dsp::basic_pvoc contract: an identity ratio reconstructs
// the input waveform exactly one frame late (the peak-locked synthesis
// guarantee), shift accuracy and level stability across the ratio range,
// clear() semantics, and float/double cross-precision agreement. The pitch
// oracle is tap::dsp::yin, certified by its own battery.

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/pvoc.h"
#include "tap/dsp/yin.h"

namespace {

    constexpr double k_pi = 3.14159265358979323846;
    constexpr double k_sr = 48000.0;

    template <typename Sample>
    std::vector<Sample> run_sine(tap::dsp::basic_pvoc<Sample>& shifter, double freq, double ratio, double seconds) {
        const int           n = static_cast<int>(seconds * k_sr);
        std::vector<Sample> out(static_cast<size_t>(n));
        for (int t = 0; t < n; ++t) {
            const Sample x              = static_cast<Sample>(std::sin(2.0 * k_pi * freq * t / k_sr));
            out[static_cast<size_t>(t)] = shifter.process(x, static_cast<Sample>(ratio));
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
    class pvoc_test : public ::testing::Test {};

    using sample_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(pvoc_test, sample_types);

    TYPED_TEST(pvoc_test, Geometry) {
        tap::dsp::basic_pvoc<TypeParam> shifter(1024);
        EXPECT_EQ(shifter.fft_size(), 1024u);
        EXPECT_EQ(shifter.hop(), 256u);
        EXPECT_EQ(shifter.latency(), 1024u);
    }

    TYPED_TEST(pvoc_test, IdentityRatioReconstructsTheWaveform) {
        // Phase-locked synthesis makes ratio 1 a true identity: the output is
        // the input delayed by exactly one FFT frame (DC and Nyquist excluded,
        // which a 440 Hz sine does not touch).
        tap::dsp::basic_pvoc<TypeParam> shifter(1024);
        const auto                      out = run_sine(shifter, 440.0, 1.0, 1.0);

        const double tolerance = std::is_same_v<TypeParam, double> ? 1e-8 : 2e-3;
        const size_t latency   = shifter.latency();
        double       worst     = 0.0;
        for (size_t t = out.size() - 4800; t < out.size(); ++t) {
            const double expected = std::sin(2.0 * k_pi * 440.0 * static_cast<double>(t - latency) / k_sr);
            worst                 = std::max(worst, std::abs(static_cast<double>(out[t]) - expected));
        }
        EXPECT_LT(worst, tolerance);
    }

    TYPED_TEST(pvoc_test, ShiftAccuracyAcrossRatioRange) {
        for (const double ratio : {0.5, 0.8, 1.122462, 1.5, 2.0}) {
            tap::dsp::basic_pvoc<TypeParam> shifter(1024);
            const auto                      out = run_sine(shifter, 220.0, ratio, 1.0);
            EXPECT_LT(std::abs(cents(measure_hz(out), 220.0 * ratio)), 8.0) << "ratio " << ratio;
            EXPECT_GT(tail_peak(out), 0.6) << "ratio " << ratio;
            EXPECT_LT(tail_peak(out), 1.4) << "ratio " << ratio;
        }
    }

    TYPED_TEST(pvoc_test, OutputStaysFinite) {
        tap::dsp::basic_pvoc<TypeParam> shifter(1024);
        for (int t = 0; t < 48000; ++t) {
            // harmonically dense input (the hard case for the class)
            const double x = 0.4 * std::sin(2.0 * k_pi * 110.0 * t / k_sr)
                             + 0.3 * std::sin(2.0 * k_pi * 347.0 * t / k_sr)
                             + 0.3 * std::sin(2.0 * k_pi * 991.0 * t / k_sr);
            const TypeParam y = shifter.process(static_cast<TypeParam>(x), TypeParam(1.26));
            ASSERT_TRUE(std::isfinite(static_cast<double>(y)));
            ASSERT_LT(std::abs(static_cast<double>(y)), 4.0);
        }
    }

    TYPED_TEST(pvoc_test, ClearZerosTheState) {
        tap::dsp::basic_pvoc<TypeParam> shifter(1024);
        run_sine(shifter, 440.0, 1.5, 0.25);
        shifter.clear();
        for (size_t i = 0; i < shifter.latency(); ++i) {
            EXPECT_EQ(shifter.process(TypeParam(0), TypeParam(1.5)), TypeParam(0));
        }
    }

    TEST(pvoc_cross_precision, FloatAgreesWithDoubleGoldenModel) {
        tap::dsp::pvoc   gold(1024);
        tap::dsp::pvoc32 fast(1024);
        const auto       od = run_sine(gold, 220.0, 1.5, 1.0);
        const auto       of = run_sine(fast, 220.0, 1.5, 1.0);
        const double     fd = measure_hz(od);
        const double     ff = measure_hz(of);
        EXPECT_NEAR(ff, fd, fd * 1e-3);
    }

} // namespace
