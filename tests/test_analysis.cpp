// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Smoke battery for the measurement instruments. These headers earn their
// keep in consumer suites (SampleRateTap's quality tests, RatioTap's alias
// acceptance); here we pin that the instruments themselves reach their
// documented floors on exact synthetic signals — an instrument that cannot
// measure a clean signal cleanly would silently weaken every consumer's
// quality gate.

#include <cmath>
#include <numbers>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/analysis/multitone_analysis.h"
#include "tap/dsp/analysis/sine_analysis.h"

namespace {

    namespace an = tap::dsp::analysis;

    TEST(SineAnalysis, FitRecoversExactTone) {
        const double       nu = 997.0 / 48000.0;
        std::vector<float> x(16384);
        for (std::size_t i = 0; i < x.size(); ++i) {
            x[i] = static_cast<float>(0.5 * std::sin(2.0 * std::numbers::pi * nu * static_cast<double>(i) + 0.3));
        }
        const an::sine_fit fit = an::fit_sine(x, nu);
        EXPECT_NEAR(fit.amplitude, 0.5, 1e-6);
        EXPECT_NEAR(fit.dc, 0.0, 1e-6);
        // Residual is the float storage quantization only: > 100 dB SNR.
        EXPECT_GT(an::snr_db(fit), 100.0);
    }

    TEST(SineAnalysis, TrackedFitAbsorbsSmallFrequencyOffset) {
        // The tone sits 5 ppm off the guess — the situation the tracker exists
        // for (a converter's rate estimate settling asymptotically). A rigid
        // fit at the guess books the offset as residual; the tracked fit must
        // recover the true frequency and the quantization-limited floor.
        const double       nu_true = (997.0 / 48000.0) * (1.0 + 5e-6);
        std::vector<float> x(65536);
        for (std::size_t i = 0; i < x.size(); ++i) {
            x[i] = static_cast<float>(0.5 * std::sin(2.0 * std::numbers::pi * nu_true * static_cast<double>(i)));
        }
        const an::sine_fit fit = an::fit_sine_tracked(x, 997.0 / 48000.0);
        EXPECT_NEAR(fit.freq_norm, nu_true, 1e-10);
        EXPECT_NEAR(fit.amplitude, 0.5, 1e-6);
        EXPECT_GT(an::snr_db(fit), 100.0);
    }

    TEST(MultitoneAnalysis, CombIsBoundedByPeakSum) {
        const an::tone_comb comb = an::tone_comb::pink(24, 40.0, 20000.0, 0.5);
        ASSERT_EQ(comb.freq_hz.size(), 24u);
        double amp_sum = 0.0;
        for (const double a : comb.amplitude) {
            amp_sum += a;
        }
        EXPECT_NEAR(amp_sum, 0.5, 1e-12);
        for (std::uint64_t i = 0; i < 48000; ++i) {
            ASSERT_LE(std::abs(comb.sample_at(i, 48000.0)), 0.5 + 1e-9);
        }
    }

    TEST(MultitoneAnalysis, JointFitReachesQuantizationFloorOnExactTones) {
        // The documented instrument floor claim: sequential fit-subtract
        // floors near 48 dB on exact synthetic tones; the joint solve reaches
        // the float storage quantization floor. Pin the latter.
        const an::tone_comb comb = an::tone_comb::pink(12, 40.0, 20000.0, 0.5);
        std::vector<float>  tail(32768);
        for (std::size_t i = 0; i < tail.size(); ++i) {
            tail[i] = static_cast<float>(comb.sample_at(i, 48000.0));
        }
        const double snr = an::program_weighted_snr_db(tail, comb, 48000.0, 48000.0);
        EXPECT_GT(snr, 90.0);
    }

} // namespace
