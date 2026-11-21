// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Pins the recursive design math (tap/dsp/iir_design.h, iir_tables.h):
// every tabulated design at the four spec points meets its spec on the
// 16384-point grid the family's designer uses, the orders reproduce the
// spike's (SampleRateTap's IIR-0, rational/tools/spikes/iir/README.md section
// 1), the constexpr coefficients sit within their pinned ulps of the mpmath
// leg (tests/reference/iir_mpmath_reference.h), and the committed tables
// promise what their rows say. The bits themselves are
// tests/test_iir_constexpr.cpp's.
//
// On the QEMU legs (tests/CMakeLists.txt, MAIN_FILTER) the grid sweeps are
// excluded by name, as the Kaiser and Nyquist searches are: dense double
// responses of every design, target-independent numbers that cost minutes
// of soft-float emulation. The structural pins (orders, ulps, symmetry, the
// committed rows' fields) run on every leg.

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

#include <gtest/gtest.h>

#include "reference/iir_mpmath_reference.h"
#include "tap/dsp/iir_design.h"

namespace {

    using namespace tap::dsp;
    namespace ref = tap::dsp::test::mpmath_ref;

    constexpr std::size_t k_grid = 16384;

    struct spec {
        int              stop_db;
        int              pass_num;
        int              pass_den;
        constexpr double p() const noexcept { return static_cast<double>(pass_num) / static_cast<double>(pass_den); }
    };
    constexpr std::array<spec, 4> k_specs{{{70, 1, 3}, {70, 3, 8}, {70, 19, 48}, {120, 5, 12}}};
    constexpr std::array<int, 4>  k_divisors{1, 2, 4, 8};

    double db20(double v) {
        return 20.0 * std::log10(std::max(v, 1e-300));
    }

    double ulps(double x, double hi, double lo) {
        const double ax  = std::fabs(x);
        const double ulp = std::nextafter(ax, std::numeric_limits<double>::infinity()) - ax;
        return std::fabs((x - hi) - lo) / ulp;
    }

    /// The worst stopband of a half-band over [1/2 - f_pass, 1/2] of its higher rate.
    template <typename D>
    double halfband_stopband_db(const D& d, double f_pass) {
        double worst = 0.0;
        for (std::size_t g = 0; g < k_grid; ++g) {
            const double f = 0.5 - f_pass + f_pass * static_cast<double>(g) / static_cast<double>(k_grid - 1);
            worst          = std::max(worst, std::abs(halfband_response(d, f)));
        }
        return db20(worst);
    }

    // ------------------------------------------------------------------------
    // The nonlinear half-band
    // ------------------------------------------------------------------------

    TEST(IirDesign, HalfbandOrdersReproduceTheSpikeAtEveryTabulatedPoint) {
        // IIR-0's table (README section 1): the fewest coefficients meeting
        // -(A + 1 dB), here from the closed form.
        constexpr std::array<std::array<int, 4>, 4> k_orders{{{9, 5, 5, 3}, {9, 5, 5, 3}, {9, 5, 5, 3}, {17, 9, 7, 5}}};
        for (std::size_t s = 0; s < k_specs.size(); ++s) {
            for (std::size_t d = 0; d < k_divisors.size(); ++d) {
                const auto des = design_nonlinear_halfband(k_specs[s].stop_db, k_specs[s].p(), k_divisors[d]);
                EXPECT_EQ(des.order, k_orders[s][d]) << k_specs[s].stop_db << " dB " << k_specs[s].pass_num << "/"
                                                     << k_specs[s].pass_den << " d" << k_divisors[d];
                EXPECT_EQ(des.count, static_cast<std::size_t>((des.order - 1) / 2));
                EXPECT_GE(des.stopband_db, k_specs[s].stop_db + k_iir_design_margin_db);
                // The next order down misses it.
                EXPECT_LT(halfband_attenuation_db(des.order - 2, des.f_pass), k_specs[s].stop_db + 1.0);
            }
        }
        static_assert(k_nonlinear_halfband<70, 3, 8, 1>.order == 9);
        static_assert(k_nonlinear_halfband<120, 5, 12, 1>.order == 17);
    }

    TEST(IirDesign, HalfbandCoefficientsAscendInTheUnitIntervalAndAlternatePaths) {
        const auto d = k_nonlinear_halfband<120, 5, 12, 1>;
        ASSERT_EQ(d.count, 8u);
        for (std::size_t i = 0; i < d.count; ++i) {
            EXPECT_GT(d.coefficients[i], 0.0);
            EXPECT_LT(d.coefficients[i], 1.0);
            if (i > 0) {
                EXPECT_GT(d.coefficients[i], d.coefficients[i - 1]);
            }
        }
        EXPECT_EQ(d.path_size(0), 4u);
        EXPECT_EQ(d.path_size(1), 4u);
        EXPECT_EQ(d.path_coefficient(1, 0), d.coefficients[1]);
        const auto n3 = k_nonlinear_halfband<70, 3, 8, 8>;
        EXPECT_EQ(n3.path_size(0), 1u);
        EXPECT_EQ(n3.path_size(1), 0u) << "an N3 half-band's second path is the bare delay";
    }

    // The tabulated designs as the compile-time constants a stage uses, in
    // the reference header's order (spec point, then divisor): their bits
    // are the same on every compiler (test_iir_constexpr.cpp), so the ulps
    // below are too, whatever the run-time arithmetic of the target.
    constexpr std::array<nonlinear_halfband_design, 16> k_tabulated_halfbands{
        k_nonlinear_halfband<70, 1, 3, 1>,   k_nonlinear_halfband<70, 1, 3, 2>,   k_nonlinear_halfband<70, 1, 3, 4>,
        k_nonlinear_halfband<70, 1, 3, 8>,   k_nonlinear_halfband<70, 3, 8, 1>,   k_nonlinear_halfband<70, 3, 8, 2>,
        k_nonlinear_halfband<70, 3, 8, 4>,   k_nonlinear_halfband<70, 3, 8, 8>,   k_nonlinear_halfband<70, 19, 48, 1>,
        k_nonlinear_halfband<70, 19, 48, 2>, k_nonlinear_halfband<70, 19, 48, 4>, k_nonlinear_halfband<70, 19, 48, 8>,
        k_nonlinear_halfband<120, 5, 12, 1>, k_nonlinear_halfband<120, 5, 12, 2>, k_nonlinear_halfband<120, 5, 12, 4>,
        k_nonlinear_halfband<120, 5, 12, 8>};

    TEST(IirDesign, HalfbandCoefficientsAreWithinFifteenUlpsOfMpmath) {
        // IIR-0 found the prototype's nome step costing up to 94,049 ulps on
        // the relaxed designs and its cancellation-free rewrite within 15
        // ulps; with it and the two further rewrites the header names (a =
        // (1 - x^2) / (1 + x)^2, and 1 - w^2 / k as a difference of squares)
        // every tabulated coefficient is within 13.2 ulps (printed below; the
        // worst is the largest coefficient at 120 dB, 5/12, d = 2, where
        // 1 - w^2 / k still cancels by about two orders of magnitude).
        ASSERT_EQ(std::size(ref::k_halfbands), k_tabulated_halfbands.size());
        double worst = 0.0;
        for (std::size_t n = 0; n < k_tabulated_halfbands.size(); ++n) {
            const ref::halfband_row&         r = ref::k_halfbands[n];
            const nonlinear_halfband_design& d = k_tabulated_halfbands[n];
            ASSERT_EQ(d.order, r.order);
            for (std::size_t i = 0; i < d.count; ++i) {
                const double e = ulps(d.coefficients[i], r.hi[i], r.lo[i]);
                worst          = std::max(worst, e);
                EXPECT_LE(e, 13.5) << r.stop_db << " dB " << r.pass_num << "/" << r.pass_den << " d" << r.divisor
                                   << " a" << i;
            }
        }
        std::printf("[ measured ] nonlinear half-band coefficients: %.2f ulps from mpmath at worst\n", worst);
        EXPECT_LE(worst, 15.0) << "IIR-0's amendment: within 15 ulps";
    }

    TEST(IirDesign, HalfbandMeetsItsSpecOnTheGrid) {
        for (const spec& s : k_specs) {
            for (const int div : k_divisors) {
                const auto   d  = design_nonlinear_halfband(s.stop_db, s.p(), div);
                const double sb = halfband_stopband_db(d, d.f_pass);
                EXPECT_LE(sb, -(s.stop_db + k_iir_design_margin_db)) << s.stop_db << " " << s.pass_num << " d" << div;
                // The grid's worst is the closed form's equiripple level.
                EXPECT_NEAR(sb, -d.stopband_db, 0.02);
            }
        }
    }

    TEST(IirDesign, HalfbandIsPowerComplementary) {
        // |H(f)|^2 + |H(1/2 - f)|^2 = 1: the bireciprocal elliptic half-band.
        for (const spec& s : k_specs) {
            const auto d     = design_nonlinear_halfband(s.stop_db, s.p(), 1);
            double     worst = 0.0;
            for (int g = 0; g <= 256; ++g) {
                const double f = 0.5 * g / 256.0;
                const double a = std::norm(halfband_response(d, f));
                const double b = std::norm(halfband_response(d, 0.5 - f));
                worst          = std::max(worst, std::fabs(a + b - 1.0));
            }
            EXPECT_LE(worst, 1e-13);
        }
    }

    TEST(IirDesign, HalfbandGroupDelayMatchesTheSpike) {
        // IIR-0, README section 2: high-rate samples, DC to the passband edge.
        struct row {
            int    stop_db, pass_num, pass_den, divisor;
            double dc, edge;
        };
        constexpr std::array<row, 4> k_rows{{{70, 1, 3, 1, 2.53, 4.21},
                                             {70, 3, 8, 1, 2.40, 4.96},
                                             {70, 19, 48, 1, 2.31, 5.54},
                                             {120, 5, 12, 1, 4.16, 11.92}}};
        for (const row& r : k_rows) {
            const auto d =
                design_nonlinear_halfband(r.stop_db, static_cast<double>(r.pass_num) / r.pass_den, r.divisor);
            EXPECT_NEAR(halfband_group_delay(d, 0.0), r.dc, 0.005);
            EXPECT_NEAR(halfband_group_delay(d, d.f_pass), r.edge, 0.005);
        }
    }

    // ------------------------------------------------------------------------
    // The Chebyshev II
    // ------------------------------------------------------------------------

    // The committed Chebyshev II designs as compile-time constants, in the
    // reference header's order (spec point, then order 3 .. 11).
    constexpr std::array<cheby2_design, 20> k_tabulated_cheby2{
        k_guard_band_cheby2<70, 1, 3, 3>,   k_guard_band_cheby2<70, 1, 3, 5>,   k_guard_band_cheby2<70, 1, 3, 7>,
        k_guard_band_cheby2<70, 1, 3, 9>,   k_guard_band_cheby2<70, 1, 3, 11>,  k_guard_band_cheby2<70, 3, 8, 3>,
        k_guard_band_cheby2<70, 3, 8, 5>,   k_guard_band_cheby2<70, 3, 8, 7>,   k_guard_band_cheby2<70, 3, 8, 9>,
        k_guard_band_cheby2<70, 3, 8, 11>,  k_guard_band_cheby2<70, 19, 48, 3>, k_guard_band_cheby2<70, 19, 48, 5>,
        k_guard_band_cheby2<70, 19, 48, 7>, k_guard_band_cheby2<70, 19, 48, 9>, k_guard_band_cheby2<70, 19, 48, 11>,
        k_guard_band_cheby2<120, 5, 12, 3>, k_guard_band_cheby2<120, 5, 12, 5>, k_guard_band_cheby2<120, 5, 12, 7>,
        k_guard_band_cheby2<120, 5, 12, 9>, k_guard_band_cheby2<120, 5, 12, 11>};

    TEST(IirDesign, Cheby2SectionsAreWithinTheirPinnedUlpsOfMpmath) {
        // Unity DC gain per section is the committed L-inf distribution at
        // every one of these (Cheby2LinfGainsAreUnityAtEveryCommittedPoint),
        // so the reference's unity-gain sections are the committed design's.
        ASSERT_EQ(std::size(ref::k_cheby2), k_tabulated_cheby2.size());
        double worst = 0.0;
        double edge  = 0.0;
        for (std::size_t n = 0; n < k_tabulated_cheby2.size(); ++n) {
            const ref::cheby2_row& r = ref::k_cheby2[n];
            const cheby2_design&   d = k_tabulated_cheby2[n];
            ASSERT_EQ(d.count, static_cast<std::size_t>((r.order + 1) / 2));
            for (std::size_t i = 0; i < d.count; ++i) {
                const iir_section&          s = d.sections[i];
                const std::array<double, 5> v{s.b0, s.b1, s.b2, s.a1, s.a2};
                for (std::size_t j = 0; j < 5; ++j) {
                    if (r.hi[5 * i + j] == 0.0) {
                        EXPECT_EQ(v[j], 0.0);
                        continue;
                    }
                    const double e = ulps(v[j], r.hi[5 * i + j], r.lo[5 * i + j]);
                    worst          = std::max(worst, e);
                    EXPECT_LE(e, 40.0) << r.stop_db << " " << r.pass_num << " order " << r.order << " s" << i << " c"
                                       << j;
                }
            }
            edge = std::max(edge, ulps(d.stop_edge, r.hi[5 * d.count], r.lo[5 * d.count]));
        }
        std::printf("[ measured ] Chebyshev II sections: %.2f ulps from mpmath at worst; stop edges %.2f\n", worst,
                    edge);
        EXPECT_LE(edge, 8.0);
    }

    TEST(IirDesign, Cheby2MeetsItsSpecWithItsZeroAtNyquist) {
        for (const spec& s : k_specs) {
            for (const int order : {3, 5, 7, 9, 11}) {
                const auto d = design_cheby2_lowpass(order, s.p(), s.stop_db);
                ASSERT_EQ(d.order, order);
                // Monotone passband: within 0.01 dB at the edge, unity at DC.
                EXPECT_NEAR(std::abs(cheby2_response(d, 0.0)), 1.0, 1e-13);
                EXPECT_GE(db20(std::abs(cheby2_response(d, s.p()))), -0.01 - 1e-9);
                // A dB from the stop edge to Nyquist, and the zero at z = -1.
                double worst = 0.0;
                for (std::size_t g = 0; g < 4096; ++g) {
                    const double f = d.stop_edge + (0.5 - d.stop_edge) * static_cast<double>(g) / 4095.0;
                    worst          = std::max(worst, std::abs(cheby2_response(d, f)));
                }
                EXPECT_LE(db20(worst), -s.stop_db + 1e-6);
                EXPECT_LE(std::abs(cheby2_response(d, 0.5)), 1e-12);
                // Ordered by increasing pole radius, every pole inside.
                for (std::size_t i = 0; i < d.count; ++i) {
                    EXPECT_LT(d.pole_radius[i], 1.0);
                    if (i > 0) {
                        EXPECT_GE(d.pole_radius[i], d.pole_radius[i - 1]);
                    }
                }
            }
        }
        // IIR-0's stop edges (README section 7), as fractions of the low Nyquist.
        EXPECT_NEAR(2.0 * design_cheby2_lowpass(5, 3.0 / 8.0, 70).stop_edge, 0.9506, 5e-5);
        EXPECT_NEAR(2.0 * design_cheby2_lowpass(7, 3.0 / 8.0, 70).stop_edge, 0.9061, 5e-5);
        EXPECT_NEAR(2.0 * design_cheby2_lowpass(5, 1.0 / 3.0, 70).stop_edge, 0.9313, 5e-5);
        EXPECT_NEAR(2.0 * design_cheby2_lowpass(7, 5.0 / 12.0, 120).stop_edge, 0.9724, 5e-5);
    }

    TEST(IirDesign, Cheby2CommittedGainsAreTheLinfGains) {
        // The committed table (tools/iir/gen_iir_tables.py, 40 digits) against
        // the same rule at run time in double: each partial cascade, ordered by
        // increasing pole radius, scaled to peak at 0 dB.
        for (const cheby2_gains& g : k_cheby2_gain_tables) {
            const auto d = design_cheby2_lowpass(g.order, static_cast<double>(g.pass_num) / g.pass_den, g.stop_db);
            std::array<double, k_max_cheby2_sections> runtime{};
            cheby2_linf_gains(d, runtime);
            const auto scaled = apply_section_gains(d, std::span<const double>(g.gains.data(), d.count - 1));
            for (std::size_t i = 0; i + 1 < d.count; ++i) {
                EXPECT_NEAR(runtime[i], g.gains[i], 1e-12);
                double peak = 0.0;
                for (int f = 0; f <= 2048; ++f) {
                    peak =
                        std::max(peak, std::abs(sections_response(
                                           std::span<const iir_section>(scaled.sections.data(), i + 1), f / 4096.0)));
                }
                EXPECT_LE(peak, 1.0 + 1e-12);
            }
            EXPECT_NEAR(std::abs(cheby2_response(scaled, 0.0)), 1.0, 1e-13);
        }
    }

    TEST(IirDesign, Cheby2LinfGainsAreUnityAtEveryCommittedPoint) {
        // A finding of this PR, pinned: with the sections ordered by
        // increasing pole radius, each at unity DC gain, every partial cascade
        // of every committed design peaks at DC, so the L-inf distribution IS
        // unity DC gain per section there (IIR-0 spent 414,830 constexpr steps
        // finding that at order 11). A custom design need not share it, which
        // is why the rule, the table and cheby2_linf_gains() stay.
        ASSERT_EQ(std::size(k_cheby2_gain_tables), 20u);
        for (const cheby2_gains& g : k_cheby2_gain_tables) {
            for (std::size_t i = 0; i < static_cast<std::size_t>((g.order - 1) / 2); ++i) {
                EXPECT_EQ(g.gains[i], 1.0);
            }
        }
        static_assert(k_guard_band_cheby2<70, 3, 8, 7>.count == 4);
        static_assert(k_guard_band_cheby2<70, 3, 8, 7>.sections[0].b0
                      == design_cheby2_lowpass(7, 3.0 / 8.0, 70).sections[0].b0);
    }

    // ------------------------------------------------------------------------
    // The committed near-linear, four-path and type II tables
    // ------------------------------------------------------------------------

    TEST(IirDesign, NearLinearTablesHoldTheSpikeOrders) {
        // IIR-0's table (README section 1): N, and K = N everywhere.
        constexpr std::array<std::array<int, 4>, 4> k_orders{{{5, 2, 2, 1}, {7, 3, 2, 1}, {8, 3, 2, 1}, {19, 5, 3, 2}}};
        ASSERT_EQ(std::size(k_near_linear_halfband_tables), 16u);
        for (std::size_t s = 0; s < k_specs.size(); ++s) {
            for (std::size_t d = 0; d < k_divisors.size(); ++d) {
                const auto* row =
                    detail::find_spec(std::span<const near_linear_halfband_design>(k_near_linear_halfband_tables),
                                      k_specs[s].stop_db, k_specs[s].pass_num, k_specs[s].pass_den, k_divisors[d]);
                ASSERT_NE(row, nullptr);
                EXPECT_EQ(row->order, k_orders[s][d]);
                EXPECT_EQ(row->delay, row->order);
                EXPECT_EQ(row->sections.size(), static_cast<std::size_t>(row->order));
                EXPECT_LT(row->pole_radius, 1.0);
                if (row->order % 2 == 1) {
                    EXPECT_GT(row->sections[row->sections.size() - 1], 0.0)
                        << "the first-order section's coefficient is positive: floor truncation leaves no limit cycle";
                }
                EXPECT_LE(row->stopband_db, -(k_specs[s].stop_db + k_iir_design_margin_db));
            }
        }
        static_assert(k_near_linear_halfband<70, 3, 8, 1>.order == 7);
        static_assert(k_near_linear_halfband<120, 5, 12, 1>.order == 19);
    }

    TEST(IirDesign, NearLinearTablesMeetTheirSpecOnTheGrid) {
        for (const near_linear_halfband_design& r : k_near_linear_halfband_tables) {
            const double f_pass = static_cast<double>(r.pass_num) / r.pass_den / (2.0 * r.divisor);
            const double sb     = halfband_stopband_db(r, f_pass);
            EXPECT_LE(sb, -(r.stop_db + k_iir_design_margin_db));
            EXPECT_NEAR(sb, r.stopband_db, 0.01);
            double dev = 0.0;
            for (std::size_t g = 0; g < k_grid; g += 4) {
                const double f = f_pass * static_cast<double>(g) / static_cast<double>(k_grid - 1);
                dev            = std::max(dev, std::fabs(halfband_group_delay(r, f) - 2.0 * r.delay));
            }
            EXPECT_LE(dev, r.delay_deviation + 1e-6);
        }
    }

    TEST(IirDesign, NearLinearPhaseErrorEqualsTheMirroredStopband) {
        // IIR-0's identity (README section 2): in the passband H = 1/2
        // e^{-j 2K w} (1 + e^{j eps}), so the phase error against the delay
        // path is eps / 2 and the stopband at the mirror is |sin(eps / 2)|.
        const auto&  r      = k_near_linear_halfband<70, 3, 8, 1>;
        const double f_pass = 3.0 / 16.0;
        double       worst  = 0.0;
        for (int g = 0; g <= 512; ++g) {
            const double f = f_pass * g / 512.0;
            const double phase =
                std::arg(halfband_response(r, f) * std::polar(1.0, 4.0 * std::numbers::pi * f * r.delay));
            const double mirror = std::abs(halfband_response(r, 0.5 - f));
            worst               = std::max(worst, std::fabs(std::fabs(phase) - std::asin(mirror)));
        }
        EXPECT_LE(worst, 1e-9);
    }

    TEST(IirDesign, FourpathTablesMeetTheFoldingBandSpec) {
        // The placements I12 admits: design divisor 2, the three 70 dB points,
        // six sections laid out [2, 2, 1, 1].
        ASSERT_EQ(std::size(k_fourpath_tables), 3u);
        for (const fourpath_design& r : k_fourpath_tables) {
            EXPECT_EQ(r.divisor, 2);
            std::size_t total = 0;
            for (const auto& p : r.paths) {
                total += p.size();
                for (const double a : p) {
                    EXPECT_GT(a, 0.0);
                    EXPECT_LT(a, 1.0);
                }
            }
            EXPECT_EQ(total, 6u);
            const double fp    = static_cast<double>(r.pass_num) / r.pass_den / (4.0 * r.divisor);
            double       worst = 0.0;
            for (const auto& band :
                 {std::array<double, 2>{0.25 - fp, 0.25 + fp}, std::array<double, 2>{0.5 - fp, 0.5}}) {
                for (std::size_t g = 0; g < k_grid / 2; ++g) {
                    const double f = band[0] + (band[1] - band[0]) * static_cast<double>(g) / (k_grid / 2 - 1);
                    worst          = std::max(worst, std::abs(fourpath_response(r, f)));
                }
            }
            EXPECT_LE(db20(worst), -(r.stop_db + k_iir_design_margin_db));
            EXPECT_NEAR(db20(worst), r.stopband_db, 0.01);
            EXPECT_NEAR(std::abs(fourpath_response(r, 0.0)), 1.0, 1e-12);
        }
        static_assert(k_fourpath<70, 3, 8, 2>.paths[0].size() == 2);
    }

    TEST(IirDesign, GuardBandFirTablesMeetTheirSpec) {
        // IIR-0's lengths (README section 7); 120 dB at order 3 did not
        // converge there either.
        ASSERT_EQ(std::size(k_guard_band_fir_tables), 19u);
        for (const guard_band_fir_design& r : k_guard_band_fir_tables) {
            const std::size_t n = r.taps.size();
            EXPECT_EQ(n % 2, 0u) << "type II: even length";
            // H(1/2) = sum_i (-1)^i h[i] pairs h[i] with h[N - 1 - i], whose
            // signs are opposite when N is even: zero, exactly.
            double alt = 0.0;
            for (std::size_t i = 0; i < n / 2; ++i) {
                EXPECT_EQ(r.taps[i], r.taps[n - 1 - i]) << "symmetric, bit for bit";
                alt += (i % 2 == 0 ? 1.0 : -1.0) * (r.taps[i] - r.taps[n - 1 - i]);
            }
            EXPECT_EQ(alt, 0.0) << "its zero at Nyquist, exactly, from its symmetry";
            const double p = static_cast<double>(r.pass_num) / r.pass_den;
            const auto   m = measure_fir_bands(r.taps, p, r.stop_edge);
            EXPECT_LE(m.ripple_db, 0.01);
            EXPECT_LE(m.stopband_db, -r.stop_db + 0.05);
            EXPECT_NEAR(r.stop_edge, design_cheby2_lowpass(r.order, p, r.stop_db).stop_edge, 1e-12);
        }
        EXPECT_EQ((k_guard_band_fir<70, 3, 8, 5>.taps.size()), 40u);
        EXPECT_EQ((k_guard_band_fir<70, 3, 8, 7>.taps.size()), 50u);
        EXPECT_EQ((k_guard_band_fir<70, 3, 8, 9>.taps.size()), 64u);
        EXPECT_EQ((k_guard_band_fir<70, 1, 3, 5>.taps.size()), 30u);
        EXPECT_EQ((k_guard_band_fir<120, 5, 12, 7>.taps.size()), 72u);
    }

    TEST(IirDesign, GuardBandLinearCustomDesignMeetsItsSpec) {
        // The run-time custom path (I4): a Kaiser-windowed type II FIR at the
        // order-5 Chebyshev's edge, 70 dB, 3/8: longer than remez's 40 taps,
        // flat within 0.01 dB, symmetric, its zero at Nyquist.
        const double      p    = 3.0 / 8.0;
        const double      edge = design_cheby2_lowpass(5, p, 70).stop_edge;
        const std::size_t n    = search_guard_band_linear_length(p, edge, 70);
        std::printf("[ measured ] Kaiser type II guard band filter, 70 dB from %.4f: %zu taps (remez: 40)\n", edge, n);
        ASSERT_GT(n, 40u);
        std::vector<double> h(n);
        design_guard_band_linear(h, p, edge, 70);
        const auto m = measure_fir_bands(h, p, edge);
        EXPECT_LE(m.stopband_db, -71.0);
        EXPECT_LE(m.ripple_db, 0.01);
        double alt = 0.0;
        for (std::size_t i = 0; i < n / 2; ++i) {
            EXPECT_EQ(h[i], h[n - 1 - i]);
            alt += (i % 2 == 0 ? 1.0 : -1.0) * (h[i] - h[n - 1 - i]);
        }
        EXPECT_EQ(alt, 0.0);
    }

    TEST(IirDesign, FloatingFlushFramesMatchTheSpike) {
        // IIR-0 (README section 5): ceil(ln(max / min_normal) / ln(1 / r)) over
        // the largest pole radius, here plus one frame per section (I6's margin).
        const double economy     = k_nonlinear_halfband<70, 3, 8, 1>.coefficients[3];   // 0.8134
        const double transparent = k_nonlinear_halfband<120, 5, 12, 1>.coefficients[7]; // 0.9188
        EXPECT_EQ(floating_flush_frames<float>(economy, 0), 853u);
        EXPECT_EQ(floating_flush_frames<double>(economy, 0), 6867u);
        EXPECT_EQ(floating_flush_frames<float>(transparent, 0), 2080u);
        EXPECT_EQ(floating_flush_frames<double>(transparent, 4), 16752u + 4u);
    }

} // namespace
