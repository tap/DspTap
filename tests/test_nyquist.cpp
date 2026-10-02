// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Locks down the tap::dsp L-th-band designer's contract (nyquist.h): the
// length arithmetic, the exact centre tap and the exact structural zeros,
// per-branch unity DC, the Nyquist property in frequency (the sum of the L
// shifted responses is 1, hence the symmetric transition f_p + f_s = r), the
// half-band MAC claim, the spec search's minimality and its measured
// numbers for the rational engine's profile candidates, and that a design
// at a non-Nyquist length is refused untouched.

#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/nyquist.h"

namespace {

    using tap::dsp::design_nyquist;
    using tap::dsp::is_nyquist_length;
    using tap::dsp::kaiser_beta;
    using tap::dsp::nyquist_centre;
    using tap::dsp::nyquist_length;
    using tap::dsp::nyquist_nonzero_taps;
    using tap::dsp::nyquist_response_db;
    using tap::dsp::nyquist_worst_stopband_db;
    using tap::dsp::search_nyquist_m;

    std::vector<double> design(std::size_t l, std::size_t m, double atten_db) {
        std::vector<double> h(nyquist_length(l, m));
        design_nyquist(h, l, kaiser_beta(atten_db));
        return h;
    }

    // Zero-phase response at f_norm = f / F_hi, referenced to the centre tap
    // (the design is symmetric, so this is real) and normalized to 1 at DC.
    double zero_phase_response(const std::vector<double>& h, std::size_t l, double f_norm) {
        const auto c   = static_cast<double>((h.size() - 1) / 2);
        double     acc = 0.0;
        for (std::size_t i = 0; i < h.size(); ++i) {
            acc += h[i] * std::cos(2.0 * std::numbers::pi * f_norm * (static_cast<double>(i) - c));
        }
        return acc / static_cast<double>(l);
    }

    TEST(Nyquist, LengthArithmetic) {
        static_assert(nyquist_length(2, 1) == 3);
        static_assert(nyquist_length(2, 9) == 35);
        static_assert(nyquist_length(3, 9) == 53);
        static_assert(nyquist_centre(2, 9) == 17);
        static_assert(nyquist_nonzero_taps(2, 9) == 19);
        static_assert(nyquist_nonzero_taps(3, 9) == 37);
        static_assert(nyquist_nonzero_taps(6, 6) == 61);
        static_assert(is_nyquist_length(2, 35));
        static_assert(!is_nyquist_length(2, 36));
        EXPECT_TRUE(is_nyquist_length(3, 53));
        EXPECT_FALSE(is_nyquist_length(3, 52));
        EXPECT_FALSE(is_nyquist_length(3, 51));
        EXPECT_FALSE(is_nyquist_length(2, 33)); // (33 + 1) % 4 != 0
        EXPECT_TRUE(is_nyquist_length(2, 3));
        EXPECT_FALSE(is_nyquist_length(2, 1));
    }

    TEST(Nyquist, CentreIsOneAndEveryLthTapIsZeroExactly) {
        for (const std::size_t l : {2u, 3u, 4u, 6u, 8u}) {
            for (const std::size_t m : {1u, 2u, 5u, 12u}) {
                const auto        h = design(l, m, 70.0);
                const std::size_t c = nyquist_centre(l, m);
                ASSERT_EQ(h.size(), nyquist_length(l, m));
                EXPECT_EQ(h[c], 1.0) << "L=" << l << " m=" << m; // exact ==
                std::size_t nonzero = 0;
                for (std::size_t i = 0; i < h.size(); ++i) {
                    const bool structural_zero = i != c && (i % l) == (c % l);
                    if (structural_zero) {
                        EXPECT_EQ(h[i], 0.0) << "L=" << l << " m=" << m << " i=" << i; // exact ==
                    }
                    else {
                        EXPECT_NE(h[i], 0.0) << "L=" << l << " m=" << m << " i=" << i;
                        ++nonzero;
                    }
                }
                EXPECT_EQ(nonzero, nyquist_nonzero_taps(l, m)) << "L=" << l << " m=" << m;
            }
        }
    }

    TEST(Nyquist, IsSymmetric) {
        const auto h = design(3, 7, 80.0);
        for (std::size_t i = 0; i < h.size(); ++i) {
            EXPECT_EQ(h[i], h[h.size() - 1 - i]) << i; // exact: the window and the sinc are both even
        }
    }

    TEST(Nyquist, EveryBranchHasUnityDcAndTheWholeSumsToL) {
        for (const std::size_t l : {2u, 3u, 6u}) {
            const auto        h     = design(l, 9, 70.0);
            const std::size_t c     = nyquist_centre(l, 9);
            double            total = 0.0;
            for (std::size_t j = 0; j < l; ++j) {
                double sum = 0.0;
                for (std::size_t i = j; i < h.size(); i += l) {
                    sum += h[i];
                }
                if (j == c % l) {
                    EXPECT_EQ(sum, 1.0) << "L=" << l << " branch " << j; // the single centre tap, exact
                }
                else {
                    EXPECT_NEAR(sum, 1.0, 1e-15)
                        << "L=" << l << " branch " << j; // double rounding only (measured 2 ulp)
                }
                total += sum;
            }
            EXPECT_NEAR(total, static_cast<double>(l), 1e-14);
        }
    }

    // The Nyquist property in frequency (Vaidyanathan 4.6): the L zero-phase
    // responses shifted by multiples of F_hi / L sum to exactly 1 at every
    // frequency, because that sum is L times the DTFT of the decimated
    // sequence h[c + nL] == delta[n]. For L = 2 this is H(f) + H(1/2 - f) = 1:
    // the response is antisymmetric about r / 2, so the transition is
    // symmetric and f_p + f_s = r.
    TEST(Nyquist, ShiftedZeroPhaseResponsesSumToOne) {
        for (const std::size_t l : {2u, 3u, 4u, 6u}) {
            const auto h = design(l, 8, 70.0);
            for (double f = 0.0; f <= 0.5; f += 0.0137) {
                double sum = 0.0;
                for (std::size_t k = 0; k < l; ++k) {
                    sum += zero_phase_response(h, l, f - static_cast<double>(k) / static_cast<double>(l));
                }
                EXPECT_NEAR(sum, 1.0, 1e-12) << "L=" << l << " f=" << f;
            }
        }
    }

    TEST(Nyquist, HalfBandTransitionIsSymmetric) {
        // H(r/2 - d) + H(r/2 + d) == 1 (zero-phase, signed): the stopband's
        // ripple is the passband's, mirrored and negated.
        const auto   h  = design(2, 9, 70.0);
        const double fc = 0.25; // r/2 of F_hi = L r = 2r
        for (double d = 0.0; d < 0.25; d += 0.02) {
            EXPECT_NEAR(zero_phase_response(h, 2, fc - d) + zero_phase_response(h, 2, fc + d), 1.0, 1e-12) << d;
        }
        EXPECT_NEAR(nyquist_response_db(h, 2, fc), -6.0206, 1e-3); // -6 dB at the lower rate's Nyquist, by symmetry
    }

    // The rational engine's profile candidates (its PLAN.md 2.3) at the
    // half-band and third-band: the searched m, its length, and the measured
    // worst stopband. Measured 2026-10-01 by this test (economy half-band
    // m = 11, N = 43, 23 nonzero, -71.9 dB; economy third-band m = 11,
    // N = 65, 45 nonzero, -71.3 dB; transparent half-band m = 31, N = 123,
    // 63 nonzero, -121.7 dB) against the plan's harris estimates of 9 / 9 /
    // 24: the Kaiser fit needs 2-7 more taps per branch for the 1 dB
    // margin. The engine's M2 spike pins the deployed counts from this
    // designer.
    struct candidate {
        std::size_t l;
        double      passband_frac; // f_p / r
        double      atten_db;
        std::size_t expected_m;
    };

    void expect_search_minimal(const candidate& c) {
        const std::size_t m = search_nyquist_m(c.l, c.passband_frac, c.atten_db);
        ASSERT_GT(m, 0u) << "L=" << c.l;
        EXPECT_EQ(m, c.expected_m) << "L=" << c.l << " p=" << c.passband_frac << " A=" << c.atten_db;
        const auto h = design(c.l, m, c.atten_db);
        EXPECT_LE(nyquist_worst_stopband_db(h, c.l, c.passband_frac), -(c.atten_db + 1.0));
        // Minimality: one branch tap fewer misses the spec (with the margin).
        if (m > 1) {
            const auto g = design(c.l, m - 1, c.atten_db);
            EXPECT_GT(nyquist_worst_stopband_db(g, c.l, c.passband_frac), -(c.atten_db + 1.0));
        }
        // The passband is flat to the edge within the Kaiser ripple.
        for (double f = 0.0; f <= c.passband_frac / static_cast<double>(c.l); f += 0.01 / static_cast<double>(c.l)) {
            EXPECT_NEAR(nyquist_response_db(h, c.l, f), 0.0, 0.02) << "L=" << c.l << " f=" << f;
        }
        std::printf("[ measured ] L=%zu p=%.4f A=%.0f dB: m=%zu N=%zu nonzero=%zu worst=%.2f dB\n", c.l,
                    c.passband_frac, c.atten_db, m, nyquist_length(c.l, m), nyquist_nonzero_taps(c.l, m),
                    nyquist_worst_stopband_db(h, c.l, c.passband_frac));
    }

    TEST(Nyquist, SearchFindsTheMinimalHalfBand) {
        // economy candidate: 70 dB, f_p = 3/8 r (18 kHz at 48 kHz).
        expect_search_minimal({.l = 2, .passband_frac = 0.375, .atten_db = 70.0, .expected_m = 11});
    }

    TEST(Nyquist, SearchFindsTheMinimalThirdBand) {
        expect_search_minimal({.l = 3, .passband_frac = 0.375, .atten_db = 70.0, .expected_m = 11});
    }

    TEST(Nyquist, SearchFindsTheMinimalTransparentHalfBand) {
        // transparent candidate: 120 dB, f_p = 5/12 r (20 kHz at 48 kHz).
        expect_search_minimal({.l = 2, .passband_frac = 5.0 / 12.0, .atten_db = 120.0, .expected_m = 31});
    }

    // The search's grid is a parameter because a long design's sidelobes
    // are narrower than the default grid's step: the 8th-band transparent
    // candidate (120 dB, p = 5/12) passes at m = 24 (N = 383) on 1024 points
    // and fails there on 16384, where m = 25 (N = 399) is the first to pass.
    // Measured 2026-10-02 (the rational engine's M2 pins use 16384 points).
    TEST(Nyquist, SearchGridIsAParameterAndResolvesLongDesignsSidelobes) {
        EXPECT_EQ(search_nyquist_m(8, 5.0 / 12.0, 120.0), 24u);
        EXPECT_EQ(search_nyquist_m(8, 5.0 / 12.0, 120.0, 1.0, 256, 16384), 25u);
        const auto h24 = design(8, 24, 120.0);
        EXPECT_LE(nyquist_worst_stopband_db(h24, 8, 5.0 / 12.0), -121.0);        // the coarse grid's verdict
        EXPECT_GT(nyquist_worst_stopband_db(h24, 8, 5.0 / 12.0, 16384), -121.0); // the sidelobe it stepped over
        EXPECT_LT(nyquist_worst_stopband_db(h24, 8, 5.0 / 12.0, 16384), -119.0);
        const auto h25 = design(8, 25, 120.0);
        EXPECT_LE(nyquist_worst_stopband_db(h25, 8, 5.0 / 12.0, 16384), -121.0);
        EXPECT_NEAR(nyquist_worst_stopband_db(h25, 8, 5.0 / 12.0, 16384),
                    nyquist_worst_stopband_db(h25, 8, 5.0 / 12.0, 65536), 0.01); // converged
    }

    TEST(Nyquist, SearchReturnsZeroWhenNoLengthMeetsTheSpec) {
        EXPECT_EQ(search_nyquist_m(2, 0.49, 120.0, 1.0, 4), 0u);
    }

    TEST(Nyquist, HalfBandComputesAboutHalfTheMacs) {
        // N = 4m - 1 taps, 2m + 1 nonzero: a decimator by 2 dots 2m + 1 of
        // them, an interpolator's branch 1 dots 2m (branch 0 is the copy).
        const std::size_t m = 9;
        EXPECT_EQ(nyquist_length(2, m), 4 * m - 1);
        EXPECT_EQ(nyquist_nonzero_taps(2, m), 2 * m + 1);
        const auto  h       = design(2, m, 70.0);
        std::size_t branch1 = 0;
        for (std::size_t i = (nyquist_centre(2, m) + 1) % 2; i < h.size(); i += 2) {
            branch1 += h[i] != 0.0;
        }
        EXPECT_EQ(branch1, 2 * m);
    }

    TEST(Nyquist, RefusesANonNyquistLengthUntouched) {
        std::vector<double> h(36, 7.0); // 2m*2 - 1 is odd; 36 is not one
        design_nyquist(h, 2, kaiser_beta(70.0));
        for (const double v : h) {
            EXPECT_EQ(v, 7.0);
        }
    }

    TEST(Nyquist, OneTapIdentityForBandOne) {
        std::vector<double> h(1, 0.0);
        design_nyquist(h, 1, kaiser_beta(70.0));
        EXPECT_EQ(h[0], 1.0);
    }

} // namespace
