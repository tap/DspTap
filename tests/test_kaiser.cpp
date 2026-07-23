// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Contract battery for the Kaiser prototype designer, ported from
// SampleRateTap's test_kaiser.cpp. The spec structs below carry that
// library's shipping preset numbers (fast/balanced/transparent/economy) so
// the designs it depends on stay pinned at their source.

#include <bit>
#include <cmath>
#include <complex>
#include <numbers>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/kaiser.h"

namespace {

    using namespace tap::dsp;

    TEST(Kaiser, BesselI0ReferenceValues) {
        EXPECT_DOUBLE_EQ(bessel_i0(0.0), 1.0);
        EXPECT_NEAR(bessel_i0(1.0), 1.2660658777520084, 1e-12);
        EXPECT_NEAR(bessel_i0(5.0), 27.239871823604442, 1e-9);
        EXPECT_NEAR(bessel_i0(12.0), 18948.925349296309, 1e-6);
    }

    TEST(Kaiser, BetaReferenceValues) {
        EXPECT_NEAR(kaiser_beta(120.0), 0.1102 * (120.0 - 8.7), 1e-12);
        EXPECT_NEAR(kaiser_beta(40.0), 0.5842 * std::pow(19.0, 0.4) + 0.07886 * 19.0, 1e-12);
        EXPECT_DOUBLE_EQ(kaiser_beta(15.0), 0.0);
    }

    TEST(Kaiser, TapEstimateMatchesHarrisFormula) {
        // 120 dB over a 20->28 kHz transition at 48 kHz: ~47 taps per phase.
        const std::size_t taps = estimate_taps(120.0, 8000.0 / 48000.0);
        EXPECT_GE(taps, 45u);
        EXPECT_LE(taps, 49u);
    }

    // A prototype specification in SampleRateTap's filter_spec vocabulary.
    struct proto_spec {
        std::size_t num_phases;
        std::size_t taps_per_phase;
        double      passband_hz;
        double      stopband_hz;
        double      stopband_atten_db;
        bool        image_zeros;

        proto_spec scaled_to(double sample_rate_hz) const {
            constexpr double k_design_rate_hz = 48000.0;
            proto_spec       s                = *this;
            s.passband_hz *= sample_rate_hz / k_design_rate_hz;
            s.stopband_hz *= sample_rate_hz / k_design_rate_hz;
            return s;
        }
    };

    // SampleRateTap's shipping presets (filter_spec::fast/balanced/transparent/
    // economy), pinned here at the substrate so a design change is caught where
    // the code now lives.
    constexpr proto_spec k_fast{128, 32, 18000.0, 30000.0, 96.0, false};
    constexpr proto_spec k_balanced{256, 48, 20000.0, 28000.0, 120.0, true};
    constexpr proto_spec k_transparent{512, 80, 20000.0, 26000.0, 140.0, true};
    constexpr proto_spec k_economy{512, 32, 18000.0, 30000.0, 96.0, true};

    // Direct DFT magnitude of the double-precision prototype, normalized so the
    // passband sits at 0 dB. f is in Hz; the prototype rate is L * fs.
    double response_db(const std::vector<double>& h, std::size_t num_phases, double fs, double f) {
        const double         proto_rate = static_cast<double>(num_phases) * fs;
        std::complex<double> acc{0.0, 0.0};
        for (std::size_t m = 0; m < h.size(); ++m) {
            const double ang = -2.0 * std::numbers::pi * f * static_cast<double>(m) / proto_rate;
            acc += h[m] * std::polar(1.0, ang);
        }
        return 20.0 * std::log10(std::abs(acc) / static_cast<double>(num_phases));
    }

    void check_prototype_meets_spec(const proto_spec& spec, double fs) {
        const std::size_t   phases = std::bit_ceil(spec.num_phases);
        const std::size_t   n      = phases * spec.taps_per_phase;
        std::vector<double> h(n);
        const double        cutoff_norm = (spec.passband_hz + spec.stopband_hz) / fs;
        if (spec.image_zeros) {
            design_prototype_compensated(h, phases, cutoff_norm, kaiser_beta(spec.stopband_atten_db),
                                         spec.passband_hz / fs);
        }
        else {
            design_prototype(h, phases, cutoff_norm, kaiser_beta(spec.stopband_atten_db));
        }

        // Passband: flat within +/-0.01 dB up to the passband edge. For the
        // compensated designs this is the claim the droop pre-compensation
        // exists to defend (the raw rect would sag -2.64 dB at 20 kHz).
        for (double f = 0.0; f <= spec.passband_hz; f += 500.0) {
            EXPECT_NEAR(response_db(h, spec.num_phases, fs, f), 0.0, 0.01) << "passband deviation at " << f << " Hz";
        }

        // Stopband: at least the rated attenuation (1 dB grace) from the stopband
        // edge out to well past the first few images.
        for (double f = spec.stopband_hz; f <= 3.0 * fs; f += 250.0) {
            EXPECT_LT(response_db(h, spec.num_phases, fs, f), -(spec.stopband_atten_db - 1.0))
                << "stopband leakage at " << f << " Hz";
        }

        // Transmission zeros at every k*fs: exact in exact arithmetic, so demand
        // far below the rated stopband (double rounding measures ~-300 dB).
        if (spec.image_zeros) {
            for (int k = 1; k <= 3; ++k) {
                EXPECT_LT(response_db(h, spec.num_phases, fs, static_cast<double>(k) * fs), -150.0)
                    << "missing transmission zero at " << k << "*fs";
            }
        }
    }

    TEST(Kaiser, FastPrototypeMeetsSpec) {
        check_prototype_meets_spec(k_fast, 48000.0);
    }

    TEST(Kaiser, BalancedPrototypeMeetsSpec) {
        check_prototype_meets_spec(k_balanced, 48000.0);
    }

    TEST(Kaiser, TransparentPrototypeMeetsSpec) {
        check_prototype_meets_spec(k_transparent, 48000.0);
    }

    TEST(Kaiser, EconomyPrototypeMeetsSpec) {
        check_prototype_meets_spec(k_economy, 48000.0);
    }

    // The compensated presets must also hold their specs at scaled rates (the
    // 16 kHz deployment path): normalized design, same numbers.
    // The k*fs transmission zeros ARE branch-DC uniformity, stated in the
    // frequency domain: with exact zeros, every polyphase branch's coefficient
    // sum is identical (measured spread 1.8e-15 -- machine epsilon -- vs 4.7e-6
    // for the plain fast design, whose spread is its stopband leakage at fs).
    TEST(Kaiser, CompensatedBranchSumsAreUniform) {
        const proto_spec    spec   = k_balanced;
        const std::size_t   phases = std::bit_ceil(spec.num_phases);
        std::vector<double> h(phases * spec.taps_per_phase);
        design_prototype_compensated(h, phases, (spec.passband_hz + spec.stopband_hz) / 48000.0,
                                     kaiser_beta(spec.stopband_atten_db), spec.passband_hz / 48000.0);
        double lo = 1e9;
        double hi = -1e9;
        for (std::size_t p = 0; p < phases; ++p) {
            double sum = 0.0;
            for (std::size_t t = 0; t < spec.taps_per_phase; ++t) {
                sum += h[t * phases + p];
            }
            lo = std::min(lo, sum);
            hi = std::max(hi, sum);
        }
        EXPECT_LT(hi - lo, 1e-12);
        EXPECT_NEAR(lo, 1.0, 1e-9);
    }

    TEST(Kaiser, CompensatedSpecsHoldAt16k) {
        check_prototype_meets_spec(k_balanced.scaled_to(16000.0), 16000.0);
        check_prototype_meets_spec(k_economy.scaled_to(16000.0), 16000.0);
    }

    // Non-power-of-two phase counts are first-class here (RatioTap designs at
    // L = 147 and 160): the designer must meet spec without the bit_ceil
    // rounding its historical consumer applied. 70 dB / 19->22.05 kHz is
    // RatioTap's economy-profile shape for the 48->44.1 direction.
    TEST(Kaiser, RationalPhaseCountMeetsSpec) {
        const std::size_t   phases = 147;
        const double        fs = 48000.0, pass = 19000.0, stop = 22050.0, atten = 70.0;
        const std::size_t   taps = estimate_taps(atten, (stop - pass) / fs);
        std::vector<double> h(phases * taps);
        design_prototype(h, phases, (pass + stop) / fs, kaiser_beta(atten));
        for (double f = 0.0; f <= pass; f += 500.0) {
            EXPECT_NEAR(response_db(h, phases, fs, f), 0.0, 0.05) << "passband deviation at " << f << " Hz";
        }
        for (double f = stop; f <= 3.0 * fs; f += 250.0) {
            EXPECT_LT(response_db(h, phases, fs, f), -(atten - 1.0)) << "stopband leakage at " << f << " Hz";
        }
    }

    TEST(Kaiser, SolveDenseSolvesKnownSystem) {
        // 3x3 that requires pivoting (zero leading entry); the solver works in
        // place, so verify by residual against saved copies.
        const std::vector<double> m0{0.0, 2.0, 1.0, /**/ 1.0, 1.0, 1.0, /**/ 2.0, 0.0, 3.0};
        const std::vector<double> rhs0{4.0, 6.0, 5.0};
        std::vector<double>       m = m0, rhs = rhs0, x(3);
        solve_dense(m, rhs, x, 3);
        for (std::size_t r = 0; r < 3; ++r) {
            const double v = m0[3 * r] * x[0] + m0[3 * r + 1] * x[1] + m0[3 * r + 2] * x[2];
            EXPECT_NEAR(v, rhs0[r], 1e-12) << "row " << r;
        }
    }

} // namespace
