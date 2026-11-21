// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Pins the floating recursive kernels (tap/dsp/iir_kernels.h), typed over
// float and double, at every tabulated design (iir_design.h, iir_tables.h):
//   - each structure equals its golden leg, the same transfer function in
//     direct form at the higher rate in double around zero-stuffing or
//     down-sampling (tests/support/iir_support.h), double within rounding and
//     float within the I14 floor;
//   - float tracks double: error power at most -137 dB re full scale at the
//     worst design (IIR-0 measured -137.1 dB at 120 dB, 5/12, near-linear,
//     interpolating, white noise; I14);
//   - DC passes at unity, a stopband tone is A dB down;
//   - the honest limit: a non-finite input poisons the state until the caller
//     zeroes it;
//   - I9's kernel-level tests, under the host's flush-to-zero set in each
//     fixture's SetUp and restored in its TearDown: silence decays to an
//     exactly zero state within floating_flush_frames(), the kernels never
//     touch the floating-point environment (the control bits, sticky flags
//     aside, equal before and after every call), and the precondition fires
//     in a Debug build with FTZ off. The no-FTZ half (the bounded residue,
//     the environment untouched with the bits clear) is
//     tests/test_iir_subnormals.cpp's, built with TAP_DSP_ALLOW_SUBNORMALS.
// Every case runs on the four QEMU legs. Where no control bit governs the
// type (double on the M4F and M33, everything on the soft-float M4),
// SilenceDecaysToExactZero states that it does not apply and passes: the
// residue there is the subnormals battery's honest limit.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numbers>
#include <span>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "support/iir_support.h"
#include "support/signals.h"
#include "tap/dsp/iir_design.h"
#include "tap/dsp/iir_kernels.h"
#include "tap/dsp/subnormals.h"

namespace {

    using namespace tap::dsp;
    using namespace tap::dsp::test;

    constexpr std::array<nonlinear_halfband_design, 16> k_nonlinear{
        k_nonlinear_halfband<70, 1, 3, 1>,   k_nonlinear_halfband<70, 1, 3, 2>,   k_nonlinear_halfband<70, 1, 3, 4>,
        k_nonlinear_halfband<70, 1, 3, 8>,   k_nonlinear_halfband<70, 3, 8, 1>,   k_nonlinear_halfband<70, 3, 8, 2>,
        k_nonlinear_halfband<70, 3, 8, 4>,   k_nonlinear_halfband<70, 3, 8, 8>,   k_nonlinear_halfband<70, 19, 48, 1>,
        k_nonlinear_halfband<70, 19, 48, 2>, k_nonlinear_halfband<70, 19, 48, 4>, k_nonlinear_halfband<70, 19, 48, 8>,
        k_nonlinear_halfband<120, 5, 12, 1>, k_nonlinear_halfband<120, 5, 12, 2>, k_nonlinear_halfband<120, 5, 12, 4>,
        k_nonlinear_halfband<120, 5, 12, 8>};

    constexpr std::array<cheby2_design, 8> k_guards{
        k_guard_band_cheby2<70, 1, 3, 5>,   k_guard_band_cheby2<70, 3, 8, 7>,   k_guard_band_cheby2<70, 19, 48, 7>,
        k_guard_band_cheby2<120, 5, 12, 7>, k_guard_band_cheby2<70, 3, 8, 11>,  k_guard_band_cheby2<70, 1, 3, 3>,
        k_guard_band_cheby2<70, 19, 48, 9>, k_guard_band_cheby2<120, 5, 12, 11>};

    constexpr std::size_t k_frames = 4096; // high-rate samples per golden comparison

    template <typename T, typename U>
    std::vector<T> cast(const std::vector<U>& x) {
        std::vector<T> y(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            y[i] = static_cast<T>(x[i]);
        }
        return y;
    }

    /// Error power of a against b, re a full-scale sine's (1/2), in dB.
    template <typename T>
    double error_db(const std::vector<T>& a, const std::vector<double>& b) {
        double e = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            const double d = static_cast<double>(a[i]) - b[i];
            e += d * d;
        }
        return 10.0 * std::log10(std::max(e / static_cast<double>(a.size()), 1e-300) / 0.5);
    }

    std::vector<double> noise(std::size_t n, std::uint32_t seed, double amplitude = 0.9) {
        xorshift32          rng(seed);
        std::vector<double> x(n);
        for (double& v : x) {
            v = amplitude * rng.next_unit();
        }
        return x;
    }

    std::vector<double> every(const std::vector<double>& x, std::size_t step) {
        std::vector<double> y;
        for (std::size_t i = 0; i < x.size(); i += step) {
            y.push_back(x[i]);
        }
        return y;
    }

    /// The pinned golden-leg error per format: double within rounding (the
    /// structures sum their rounding differently; measured -309.4 dB at worst
    /// on x86-64), float the I14 floor (measured -138.1 dB at worst, the
    /// near-linear interpolators, against IIR-0's -137.1).
    template <typename T>
    constexpr double k_golden_db = std::is_same_v<T, double> ? -300.0 : -137.0;

    template <typename T>
    class iir_kernel : public flush_to_zero_fixture {};

    using floating_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(iir_kernel, floating_types);

    TYPED_TEST(iir_kernel, HalfbandDecimatorMatchesTheDirectFormGoldenLeg) {
        using real   = TypeParam;
        double worst = -400.0;
        for (const auto& d : k_nonlinear) {
            const std::vector<double> x   = cast<double>(cast<real>(noise(k_frames, 0x1234u)));
            const std::vector<double> ref = every(halfband_direct(d, x), 2);
            const auto                y   = run_halfband_decimate(nonlinear_halfband_table<real>(d), cast<real>(x));
            const double              e   = error_db(y, ref);
            worst                         = std::max(worst, e);
            EXPECT_LE(e, k_golden_db<real>) << "order " << d.order << " f_pass " << d.f_pass;
        }
        std::printf("[ measured ] nonlinear decimator vs the direct form: %.1f dB at worst\n", worst);
    }

    TYPED_TEST(iir_kernel, HalfbandInterpolatorMatchesTheDirectFormGoldenLeg) {
        using real   = TypeParam;
        double worst = -400.0;
        for (const auto& d : k_nonlinear) {
            const std::vector<double> x   = cast<double>(cast<real>(noise(k_frames / 2, 0x2345u)));
            std::vector<double>       ref = halfband_direct(d, zero_stuffed(x, 2));
            for (double& v : ref) {
                v *= 2.0;
            }
            const auto   y = run_halfband_interpolate(nonlinear_halfband_table<real>(d), cast<real>(x));
            const double e = error_db(y, ref);
            worst          = std::max(worst, e);
            EXPECT_LE(e, k_golden_db<real>) << "order " << d.order;
        }
        std::printf("[ measured ] nonlinear interpolator vs the direct form: %.1f dB at worst\n", worst);
    }

    TYPED_TEST(iir_kernel, NearLinearKernelsMatchTheDirectFormGoldenLeg) {
        using real     = TypeParam;
        double worst_d = -400.0;
        double worst_i = -400.0;
        for (const auto& d : k_near_linear_halfband_tables) {
            const near_linear_halfband_table<real> t(d);
            const std::vector<double>              x   = cast<double>(cast<real>(noise(k_frames, 0x3456u)));
            const std::vector<double>              ref = every(halfband_direct(d, x), 2);
            worst_d                      = std::max(worst_d, error_db(run_near_linear_decimate(t, cast<real>(x)), ref));
            const std::vector<double> xi = cast<double>(cast<real>(noise(k_frames / 2, 0x4567u)));
            std::vector<double>       refi = halfband_direct(d, zero_stuffed(xi, 2));
            for (double& v : refi) {
                v *= 2.0;
            }
            worst_i = std::max(worst_i, error_db(run_near_linear_interpolate(t, cast<real>(xi)), refi));
        }
        std::printf("[ measured ] near-linear decimator / interpolator vs the direct form: %.1f / %.1f dB\n", worst_d,
                    worst_i);
        EXPECT_LE(worst_d, k_golden_db<real>);
        EXPECT_LE(worst_i, k_golden_db<real>);
    }

    TYPED_TEST(iir_kernel, FourpathKernelsMatchTheDirectFormGoldenLeg) {
        using real     = TypeParam;
        double worst_d = -400.0;
        double worst_i = -400.0;
        for (const auto& d : k_fourpath_tables) {
            const fourpath_table<real> t(d);
            const std::vector<double>  x   = cast<double>(cast<real>(noise(k_frames, 0x5678u)));
            const std::vector<double>  ref = every(fourpath_direct(d, x), 4);
            worst_d                        = std::max(worst_d, error_db(run_fourpath_decimate(t, cast<real>(x)), ref));
            const std::vector<double> xi   = cast<double>(cast<real>(noise(k_frames / 4, 0x6789u)));
            std::vector<double>       refi = fourpath_direct(d, zero_stuffed(xi, 4));
            for (double& v : refi) {
                v *= 4.0;
            }
            worst_i = std::max(worst_i, error_db(run_fourpath_interpolate(t, cast<real>(xi)), refi));
        }
        std::printf("[ measured ] four-path decimator / interpolator vs the direct form: %.1f / %.1f dB\n", worst_d,
                    worst_i);
        EXPECT_LE(worst_d, k_golden_db<real>);
        EXPECT_LE(worst_i, k_golden_db<real>);
    }

    TYPED_TEST(iir_kernel, BiquadCascadeMatchesTheDirectFormGoldenLeg) {
        // Transposed direct form II against direct form I: the same transfer
        // function, a different rounding.
        using real   = TypeParam;
        double worst = -400.0;
        for (const auto& d : k_guards) {
            const std::vector<double> x = cast<double>(cast<real>(noise(k_frames, 0x789au)));
            const std::vector<double> ref =
                sections_direct(std::span<const iir_section>(d.sections.data(), d.count), x);
            worst = std::max(worst, error_db(run_biquad_cascade(biquad_cascade_table<real>(d), cast<real>(x)), ref));
        }
        // Measured on x86-64: double -296.5 dB; float -125.4 dB at worst (120
        // dB, 5/12, order 11, its poles out to radius 0.952), -126.1 at 120
        // dB order 7, -130.3 to -136.7 dB at the 70 dB points: float states
        // and products (I14) hold the guard band filter above its spec, by 5
        // dB at 120 dB. The pin leaves a dB for targets that contract.
        std::printf("[ measured ] Chebyshev II cascade vs direct form I: %.1f dB at worst\n", worst);
        constexpr double pin = std::is_same_v<real, double> ? -290.0 : -124.0;
        EXPECT_LE(worst, pin);
    }

    TYPED_TEST(iir_kernel, DcPassesAtUnity) {
        using real          = TypeParam;
        const auto&       d = k_nonlinear[4];
        const std::size_t n = 4000;
        const auto        y = run_halfband_decimate(nonlinear_halfband_table<real>(d), std::vector<real>(n, real(0.5)));
        const auto   yi  = run_halfband_interpolate(nonlinear_halfband_table<real>(d), std::vector<real>(n, real(0.5)));
        const auto   g   = run_biquad_cascade(biquad_cascade_table<real>(k_guards[1]), std::vector<real>(n, real(0.5)));
        const double tol = std::is_same_v<real, double> ? 1e-14 : 1e-6;
        EXPECT_NEAR(static_cast<double>(y.back()), 0.5, tol);
        EXPECT_NEAR(static_cast<double>(yi.back()), 0.5, tol);
        EXPECT_NEAR(static_cast<double>(yi[yi.size() - 2]), 0.5, tol);
        EXPECT_NEAR(static_cast<double>(g.back()), 0.5, tol);
    }

    TYPED_TEST(iir_kernel, AStopbandToneIsAttenuatedByTheSpec) {
        // A tone at 31/64 of the higher rate, in every tabulated stopband,
        // decimated: its alias at 31/32 of the output rate has, in steady
        // state, the amplitude |H(31/64)| / 2 the design's response states,
        // at least A + 1 dB down (float: down to its own floor at 120 dB).
        // The amplitude is the RMS over whole periods (32 outputs) times
        // sqrt 2, exact for a sampled sinusoid.
        using real            = TypeParam;
        constexpr double    f = 31.0 / 64.0;
        std::vector<double> x(16384);
        for (std::size_t i = 0; i < x.size(); ++i) {
            x[i] = 0.5 * std::cos(2.0 * std::numbers::pi * f * static_cast<double>(i));
        }
        for (const auto& d : k_nonlinear) {
            ASSERT_GE(f, 0.5 - d.f_pass);
            const auto y   = run_halfband_decimate(nonlinear_halfband_table<real>(d), cast<real>(x));
            double     pow = 0.0;
            for (std::size_t i = y.size() / 2; i < y.size(); ++i) {
                pow += static_cast<double>(y[i]) * static_cast<double>(y[i]);
            }
            const double amplitude = std::sqrt(2.0 * pow / static_cast<double>(y.size() / 2));
            // The float output's own rounding floor bounds what float can show.
            constexpr double floor  = std::is_same_v<real, double> ? 1e-15 : 3e-8;
            const double     expect = 0.5 * std::abs(halfband_response(d, f));
            EXPECT_NEAR(amplitude, expect, floor + 1e-9 * expect) << "order " << d.order;
            EXPECT_LE(amplitude, 0.5 * std::pow(10.0, -(d.stopband_db - 0.01) / 20.0) + floor);
        }
    }

    TYPED_TEST(iir_kernel, NonFiniteInputPoisonsTheStateUntilItIsZeroed) {
        // The honest limit (I9, decided): an FIR clears a non-finite burst
        // after one window, a recursion never does; zeroing the state is the
        // way out, and after it the kernel is a fresh one, bit for bit.
        using real = TypeParam;
        const nonlinear_halfband_table<real> t(k_nonlinear[4]);
        std::vector<iir_state<real>>         st(t.state_size(), iir_state<real>{});
        const std::vector<real>              x = cast<real>(noise(512, 0xabcdu));
        static_cast<void>(halfband_decimate(t, st.data(), std::numeric_limits<real>::quiet_NaN(), real(0)));
        for (std::size_t k = 0; k < 256; ++k) {
            EXPECT_TRUE(std::isnan(halfband_decimate(t, st.data(), x[2 * k], x[2 * k + 1]))) << "frame " << k;
        }
        std::fill(st.begin(), st.end(), iir_state<real>{});
        std::vector<iir_state<real>> fresh(t.state_size(), iir_state<real>{});
        for (std::size_t k = 0; k < 256; ++k) {
            EXPECT_EQ(halfband_decimate(t, st.data(), x[2 * k], x[2 * k + 1]),
                      halfband_decimate(t, fresh.data(), x[2 * k], x[2 * k + 1]));
        }
        // An infinity poisons it the same way (inf - inf is NaN inside).
        static_cast<void>(halfband_decimate(t, st.data(), std::numeric_limits<real>::infinity(), real(0)));
        bool finite_again = false;
        for (std::size_t k = 0; k < 256; ++k) {
            finite_again = finite_again || std::isfinite(halfband_decimate(t, st.data(), real(0), real(0)));
        }
        EXPECT_FALSE(finite_again);
    }

    // ------------------------------------------------------------------------
    // I9, the kernels' counterparts of the plan's IirStage tests
    // ------------------------------------------------------------------------

    /// The largest pole radius per path sample of a design: its largest
    /// one-multiplier coefficient.
    [[maybe_unused]] double largest_radius(const nonlinear_halfband_design& d) {
        return d.coefficients[d.count - 1];
    }

    TYPED_TEST(iir_kernel, SilenceDecaysToExactZero) {
        using real = TypeParam;
        if constexpr (!k_subnormal_control<real>) {
            // No bit governs real here (software float): nothing flushes, and the
            // residue is test_iir_subnormals.cpp's honest limit.
            EXPECT_FALSE(subnormals_flushed<real>());
            return;
        }
        else {
            ASSERT_TRUE(subnormals_flushed<real>()) << "the fixture set the host's flush-to-zero";
            const auto zero = [](const std::vector<iir_state<real>>& st) {
                for (const auto s : st) {
                    if (s != iir_state<real>{}) {
                        return false;
                    }
                }
                return true;
            };
            const std::vector<real> x       = cast<real>(noise(4000, 0x1111u, 1.0));
            std::size_t             longest = 0;
            for (const auto& d : k_nonlinear) {
                const nonlinear_halfband_table<real> t(d);
                const std::size_t            window = floating_flush_frames<real>(largest_radius(d), t.multiplies());
                std::vector<iir_state<real>> st(t.state_size(), iir_state<real>{});
                real                         odd{};
                for (std::size_t k = 0; k < 2000; ++k) {
                    static_cast<void>(halfband_decimate(t, st.data(), x[2 * k], odd));
                    odd = x[2 * k + 1];
                }
                std::size_t frames = 0;
                for (; !zero(st) && frames < 4 * window; ++frames) {
                    static_cast<void>(halfband_decimate(t, st.data(), real(0), odd));
                    odd = real(0);
                }
                EXPECT_TRUE(zero(st));
                EXPECT_LE(frames, window) << "order " << d.order << ": the state is exactly zero within the window";
                EXPECT_EQ(halfband_decimate(t, st.data(), real(0), real(0)), real(0)) << "and the output is silence";
                longest = std::max(longest, frames);
            }
            for (const auto& d : k_near_linear_halfband_tables) {
                const near_linear_halfband_table<real> t(d);
                const std::size_t            window = floating_flush_frames<real>(d.pole_radius, t.multiplies());
                std::vector<iir_state<real>> st(t.state_size(), iir_state<real>{});
                for (std::size_t k = 0; k < 2000; ++k) {
                    static_cast<void>(near_linear_decimate(t, st.data(), x[2 * k], x[2 * k + 1]));
                }
                std::size_t frames = 0;
                for (; !zero(st) && frames < 4 * window; ++frames) {
                    static_cast<void>(near_linear_decimate(t, st.data(), real(0), real(0)));
                }
                EXPECT_TRUE(zero(st));
                EXPECT_LE(frames, window) << "near-linear N" << d.order;
            }
            std::printf("[ measured ] %s: the nonlinear states reach exactly zero within %lu frames at worst\n",
                        std::is_same_v<real, float> ? "float" : "double", static_cast<unsigned long>(longest));
        }
    }

    TYPED_TEST(iir_kernel, TheGuardBandCascadeHoldsALimitCycleAtTheSmallestNormal) {
        // The honest limit of the floating biquad cascade (iir_kernels.h):
        // under flush-to-zero, once |a2 y| falls below the smallest normal its
        // product flushes and the section runs for a sample as a first-order
        // recursion with |a1| > 1, which lifts it again; the Chebyshev II of
        // order 7 or more holds that zero-input limit cycle at a few times the
        // smallest normal and never reaches an exactly zero state (direct form
        // I does the same: measured 3.8 to 7.1 times, transposed direct form II
        // 4.0 to 13.1). Orders 3 and 5 at 1/3 drain to zero within I6's count.
        using real = TypeParam;
        if constexpr (!k_subnormal_control<real>) {
            EXPECT_FALSE(subnormals_flushed<real>());
            return;
        }
        else {
            // Driven by a full-scale square wave of the given period, then
            // silent for `silence` samples: the largest state left, in units
            // of the smallest normal.
            const auto run = [](const cheby2_design& d, int period, std::size_t silence) {
                const biquad_cascade_table<real> t(d);
                std::vector<iir_state<real>>     st(t.state_size(), iir_state<real>{});
                for (int k = 0; k < 3000; ++k) {
                    static_cast<void>(biquad_cascade(t, st.data(), (k % period) < period / 2 ? real(1) : real(-1)));
                }
                for (std::size_t n = 0; n < silence; ++n) {
                    static_cast<void>(biquad_cascade(t, st.data(), real(0)));
                }
                double worst = 0.0;
                for (const auto s : st) {
                    worst = std::max(worst, std::fabs(static_cast<double>(s)));
                }
                return worst / static_cast<double>(std::numeric_limits<real>::min());
            };
            // Whether a start ends in the cycle depends on the state it starts
            // from; across these, some do at each design, and none exceeds the
            // bound.
            for (const auto& d : {k_guard_band_cheby2<70, 3, 8, 7>, k_guard_band_cheby2<120, 5, 12, 11>}) {
                const std::size_t silence = 20 * floating_flush_frames<real>(d.pole_radius[d.count - 1], d.count);
                double            worst   = 0.0;
                int               cycling = 0;
                for (const int period : {6, 10, 14, 22, 37}) {
                    const double residue = run(d, period, silence);
                    worst                = std::max(worst, residue);
                    cycling += residue > 0.0 ? 1 : 0;
                }
                std::printf("[ measured ] %s Chebyshev II order %d: after 20 windows of silence %d of 5 starts hold a "
                            "cycle, at most %.1f x the smallest normal\n",
                            std::is_same_v<real, float> ? "float" : "double", d.order, cycling, worst);
                EXPECT_GT(cycling, 0) << "the limit cycle the header states";
                EXPECT_LE(worst, 16.0);
            }
            const auto& d5 = k_guard_band_cheby2<70, 1, 3, 5>;
            for (const int period : {6, 10, 14, 22, 37}) {
                EXPECT_EQ(run(d5, period, floating_flush_frames<real>(d5.pole_radius[d5.count - 1], d5.count)), 0.0)
                    << "order 5 at 1/3 drains to an exactly zero state";
            }
        }
    }

    TYPED_TEST(iir_kernel, NeverTouchesTheFloatingPointEnvironment) {
        // The kernels only read the control register (the debug-build
        // precondition); its control bits are the same after every call.
        using real = TypeParam;
        const nonlinear_halfband_table<real>   hb(k_nonlinear[4]);
        const near_linear_halfband_table<real> nl(k_near_linear_halfband<70, 3, 8, 1>);
        const fourpath_table<real>             fp(k_fourpath<70, 3, 8, 2>);
        const biquad_cascade_table<real>       bq(k_guards[1]);
        std::vector<iir_state<real>> s1(hb.state_size()), s2(nl.state_size()), s3(fp.state_size()), s4(bq.state_size());
        const std::vector<real>      x = cast<real>(noise(64, 0x2222u));
        std::array<real, 4>          out{};
        for (std::size_t i = 0; i + 4 <= x.size(); i += 4) {
            const std::uint64_t before = floating_point_control() & k_control_mask;
            static_cast<void>(halfband_decimate(hb, s1.data(), x[i], x[i + 1]));
            halfband_interpolate(hb, s1.data(), x[i], out.data());
            static_cast<void>(near_linear_decimate(nl, s2.data(), x[i + 2], x[i + 1]));
            near_linear_interpolate(nl, s2.data(), x[i + 2], x[i], out.data());
            static_cast<void>(fourpath_decimate(fp, s3.data(), x.data() + i));
            fourpath_interpolate(fp, s3.data(), x[i], out.data());
            static_cast<void>(biquad_cascade(bq, s4.data(), x[i + 3]));
            EXPECT_EQ(floating_point_control() & k_control_mask, before);
        }
    }

    // ------------------------------------------------------------------------
    // Cross-precision: float states against double (I14)
    // ------------------------------------------------------------------------

    class iir_kernel_cross_precision : public flush_to_zero_fixture {};

    TEST_F(iir_kernel_cross_precision, FloatTracksDoubleWithinTheI14Floor) {
        // White noise at 0.9 of full scale through every tabulated structure,
        // float against double on the same float input; error power re a
        // full-scale sine's. IIR-0 measured -137.1 dB at worst (120 dB, 5/12,
        // near-linear, interpolating), -141.7 dB or lower at 70 dB.
        double worst = -400.0;
        double w70   = -400.0;
        for (const auto& d : k_nonlinear) {
            const std::vector<float> x  = cast<float>(noise(1 << 14, 0x3333u));
            const auto               yf = run_halfband_decimate(nonlinear_halfband_table<float>(d), x);
            const auto   yd = run_halfband_decimate(nonlinear_halfband_table<double>(d), cast<double>(cast<double>(x)));
            const double e  = error_db(yf, yd);
            const auto   yfi = run_halfband_interpolate(nonlinear_halfband_table<float>(d), x);
            const auto   ydi = run_halfband_interpolate(nonlinear_halfband_table<double>(d), cast<double>(x));
            const double ei  = error_db(yfi, ydi);
            worst            = std::max({worst, e, ei});
            if (d.order <= 9 && d.stopband_db < 100.0) {
                w70 = std::max({w70, e, ei});
            }
        }
        for (const auto& d : k_near_linear_halfband_tables) {
            const std::vector<float> x  = cast<float>(noise(1 << 14, 0x4444u));
            const auto               yf = run_near_linear_decimate(near_linear_halfband_table<float>(d), x);
            const auto yd  = run_near_linear_decimate(near_linear_halfband_table<double>(d), cast<double>(x));
            const auto yfi = run_near_linear_interpolate(near_linear_halfband_table<float>(d), x);
            const auto ydi = run_near_linear_interpolate(near_linear_halfband_table<double>(d), cast<double>(x));
            worst          = std::max({worst, error_db(yf, yd), error_db(yfi, ydi)});
            if (d.stop_db == 70) {
                w70 = std::max({w70, error_db(yf, yd), error_db(yfi, ydi)});
            }
        }
        std::printf("[ measured ] float states against double: %.1f dB at worst, %.1f dB at the 70 dB points\n", worst,
                    w70);
        EXPECT_LE(worst, -137.0);
        EXPECT_LE(w70, -140.0);
    }

#if !defined(NDEBUG) && defined(GTEST_HAS_DEATH_TEST) && GTEST_HAS_DEATH_TEST
    // The precondition itself: a debug assertion (TAP_EXPECTS), compiled in
    // only where it can fire and be caught (a Debug configure on a host with
    // death tests, as fft_engine_debug's), and only where a bit governs float.
    TEST(IirKernelDebug, ThePreconditionFiresWithoutFlushToZero) {
        if constexpr (k_subnormal_control<float>) {
            const nonlinear_halfband_table<float> t(k_nonlinear[4]);
            std::vector<iir_state<float>>         st(t.state_size());
            const std::uint64_t                   saved = set_flush_to_zero(false);
            EXPECT_DEATH(static_cast<void>(halfband_decimate(t, st.data(), 0.5f, 0.0f)), "");
            write_floating_point_control(saved);
            const std::uint64_t on = set_flush_to_zero(true);
            EXPECT_EQ(halfband_decimate(t, st.data(), 0.0f, 0.0f), 0.0f) << "with the bits set it does not fire";
            write_floating_point_control(on);
        }
    }
#endif

} // namespace
