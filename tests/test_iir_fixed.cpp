// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Pins the fixed-point recursive kernels (tap/dsp/iir_kernels.h, IIR_PLAN I8),
// typed over Q15 (int16_t) and Q31 (int32_t), at every tabulated design:
//   - unity gain: a constant passes exactly (an allpass is exactly unity at
//     DC, and the truncated recursions settle on the input), the guard band
//     cascade within one LSB (its quantized numerators);
//   - the headroom: full-scale square waves and clipped sines never take a
//     state past two integer guard bits, and every output stays within one
//     LSB of the double structure rounded (so nothing wrapped); the largest
//     state is printed in units of full scale;
//   - saturation at the output only: an overshoot past full scale clamps;
//   - the limit-cycle decay: silence after full-scale drive reaches an exactly
//     zero state within the stated count, the nonlinear structures' bound
//     computed here exactly (each one-multiplier section's longest floor-
//     truncated decay from the full state headroom, summed along each path,
//     IIR-0's method), the near-linear and guard band ones pinned as measured
//     beside IIR-0's magnitude-truncation bounds;
//   - the Q15 floor: a -1 dBFS sine through the Q15 decimator is within
//     0.1 dB of the double structure's output rounded to Q15 (IIR-0: -98.1 dB
//     at 70 dB, 3/8, the recursion adding nothing above the output's own
//     rounding);
//   - the bits: FNV-1a-64 of every structure's Q15 and Q31 outputs on a fixed
//     input, one row for every compiler and target CI runs
//     (tests/support/iir_fixed_pins.h): integer arithmetic, no libm, the
//     coefficient tables bit-pinned by tests/test_iir_constexpr.cpp.
// Every case runs on the four QEMU legs.

#include <algorithm>
#include <array>
#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numbers>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "support/iir_fixed_pins.h"
#include "support/iir_support.h"
#include "support/signals.h"
#include "tap/dsp/iir_design.h"
#include "tap/dsp/iir_kernels.h"

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

    constexpr std::array<cheby2_design, 4> k_guards{k_guard_band_cheby2<70, 1, 3, 5>, k_guard_band_cheby2<70, 3, 8, 7>,
                                                    k_guard_band_cheby2<70, 19, 48, 7>,
                                                    k_guard_band_cheby2<120, 5, 12, 7>};

    template <typename S>
    constexpr double k_scale = static_cast<double>(std::int64_t{1} << std::numeric_limits<S>::digits);

    template <typename S>
    std::vector<S> quantize(const std::vector<double>& x) {
        std::vector<S> y(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            y[i] = sample_scale<S>::from_double(x[i]);
        }
        return y;
    }

    template <typename S>
    std::vector<double> unit(const std::vector<S>& x) {
        std::vector<double> y(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            y[i] = sample_scale<S>::to_double(x[i]);
        }
        return y;
    }

    /// The largest |state| of a state vector, in units of full scale.
    template <typename S>
    double largest_state(const std::vector<iir_state<S>>& st) {
        double m = 0.0;
        for (const auto v : st) {
            m = std::max(m, std::fabs(static_cast<double>(v)));
        }
        return m / static_cast<double>(std::int64_t{1} << iir_traits<S>::k_state_frac_bits);
    }

    /// The test inputs of IIR-0's headroom check: full-scale square waves
    /// (periods 2 to 400 high-rate samples, a spread of them) and sines at
    /// +12 dB clipped to full scale.
    std::vector<std::vector<double>> full_scale_inputs(std::size_t n) {
        std::vector<std::vector<double>> out;
        for (const int per : {2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 32, 48, 64, 100, 150, 200, 300, 400}) {
            std::vector<double> x(n);
            for (std::size_t i = 0; i < n; ++i) {
                x[i] = (static_cast<int>(i % static_cast<std::size_t>(per)) < per / 2) ? 1.0 : -1.0;
            }
            out.push_back(x);
        }
        for (const double fr : {0.001, 0.01, 0.03, 0.07, 0.11, 0.17, 0.23, 0.29}) {
            std::vector<double> x(n);
            for (std::size_t i = 0; i < n; ++i) {
                x[i] = std::clamp(4.0 * std::sin(2.0 * std::numbers::pi * fr * static_cast<double>(i)), -1.0, 1.0);
            }
            out.push_back(x);
        }
        return out;
    }

    template <typename S>
    class iir_fixed : public ::testing::Test {};

    using fixed_types = ::testing::Types<std::int16_t, std::int32_t>;
    TYPED_TEST_SUITE(iir_fixed, fixed_types);

    TYPED_TEST(iir_fixed, AConstantPassesAtUnityExactly) {
        using sample = TypeParam;
        for (const double level : {1.0, -1.0, 0.5, -0.3, 1e-3}) {
            const sample              c = sample_scale<sample>::from_double(level);
            const std::vector<sample> x(4000, c);
            for (const auto& d : k_nonlinear) {
                const nonlinear_halfband_table<sample> t(d);
                EXPECT_EQ(run_halfband_decimate(t, x).back(), c) << "order " << d.order << " level " << level;
                const auto yi = run_halfband_interpolate(t, x);
                EXPECT_EQ(yi.back(), c);
                EXPECT_EQ(yi[yi.size() - 2], c);
            }
            for (const auto& d : k_near_linear_halfband_tables) {
                const near_linear_halfband_table<sample> t(d);
                EXPECT_EQ(run_near_linear_decimate(t, x).back(), c) << "near-linear N" << d.order;
                EXPECT_EQ(run_near_linear_interpolate(t, x).back(), c);
            }
            for (const auto& d : k_fourpath_tables) {
                EXPECT_EQ(run_fourpath_decimate(fourpath_table<sample>(d), x).back(), c);
                EXPECT_EQ(run_fourpath_interpolate(fourpath_table<sample>(d), x).back(), c);
            }
            for (const auto& d : k_guards) {
                const auto y = run_biquad_cascade(biquad_cascade_table<sample>(d), x);
                EXPECT_LE(std::fabs(static_cast<double>(y.back()) - static_cast<double>(c)), 1.0)
                    << "the guard band cascade within one LSB of unity (its Q2.30 numerators)";
            }
        }
    }

    // ------------------------------------------------------------------------
    // The headroom: two guard bits hold every state, nothing saturates inside
    // ------------------------------------------------------------------------

    // The double structure on the format's own quantized coefficients: the
    // reference that isolates the fixed-point arithmetic (at a Q31 output the
    // Q0.31 / Q2.30 coefficients' own quantization is visible, several LSB on
    // the guard band cascades; at Q15 it is not).
    template <typename S>
    constexpr double coeff1_value(typename iir_traits<S>::coeff1 c) {
        return static_cast<double>(c) / static_cast<double>(std::int64_t{1} << iir_traits<S>::k_coeff1_frac_bits);
    }
    template <typename S>
    constexpr double coeff2_value(typename iir_traits<S>::coeff2 c) {
        return static_cast<double>(c) / static_cast<double>(std::int64_t{1} << iir_traits<S>::k_coeff2_frac_bits);
    }
    template <typename S>
    nonlinear_halfband_table<double> as_double(const nonlinear_halfband_table<S>& t) {
        nonlinear_halfband_table<double> d;
        d.size0 = t.size0;
        d.size1 = t.size1;
        for (std::size_t i = 0; i < t.size0; ++i) {
            d.path0[i] = coeff1_value<S>(t.path0[i]);
        }
        for (std::size_t i = 0; i < t.size1; ++i) {
            d.path1[i] = coeff1_value<S>(t.path1[i]);
        }
        return d;
    }
    template <typename S>
    near_linear_halfband_table<double> as_double(const near_linear_halfband_table<S>& t) {
        near_linear_halfband_table<double> d;
        d.n2    = t.n2;
        d.n1    = t.n1;
        d.delay = t.delay;
        for (std::size_t i = 0; i < 2 * t.n2; ++i) {
            d.second[i] = coeff2_value<S>(t.second[i]);
        }
        d.first[0] = coeff1_value<S>(t.first[0]);
        return d;
    }
    template <typename S>
    fourpath_table<double> as_double(const fourpath_table<S>& t) {
        fourpath_table<double> d;
        d.sizes = t.sizes;
        for (std::size_t k = 0; k < 4; ++k) {
            for (std::size_t i = 0; i < t.sizes[k]; ++i) {
                d.paths[k][i] = coeff1_value<S>(t.paths[k][i]);
            }
        }
        return d;
    }
    template <typename S>
    biquad_cascade_table<double> as_double(const biquad_cascade_table<S>& t) {
        biquad_cascade_table<double> d;
        d.size = t.size;
        for (std::size_t i = 0; i < t.size; ++i) {
            for (std::size_t j = 0; j < 5; ++j) {
                d.sections[i][j] = coeff2_value<S>(t.sections[i][j]);
            }
        }
        return d;
    }

    // Each structure runs in the format and in double (on the same quantized
    // coefficients) side by side; the largest state (units of full scale) and
    // the largest output error (LSB, against the double output saturated to
    // the format) are tracked.
    TYPED_TEST(iir_fixed, FullScaleInputNeverSaturatesInside) {
        using sample             = TypeParam;
        const auto   inputs      = full_scale_inputs(1200);
        double       worst_state = 0.0;
        double       worst_err   = 0.0;
        const double lsb         = 1.0 / k_scale<sample>;
        const auto   track       = [&](double state, double err) {
            worst_state = std::max(worst_state, state);
            worst_err   = std::max(worst_err, err);
        };
        for (const auto& x : inputs) {
            const std::vector<sample> xq = quantize<sample>(x);
            const std::vector<double> xd = unit(xq);
            for (const auto& d : k_nonlinear) {
                const nonlinear_halfband_table<sample> tq(d);
                const nonlinear_halfband_table<double> td = as_double(tq);
                std::vector<iir_state<sample>>         sq(tq.state_size(), iir_state<sample>{});
                std::vector<double>                    sd(td.state_size(), 0.0);
                sample                                 oq{};
                double                                 od = 0.0;
                for (std::size_t k = 0; k + 1 < x.size(); k += 2) {
                    const sample yq = halfband_decimate(tq, sq.data(), xq[k], oq);
                    const double yd = halfband_decimate(td, sd.data(), xd[k], od);
                    oq              = xq[k + 1];
                    od              = xd[k + 1];
                    track(largest_state<sample>(sq),
                          std::fabs(sample_scale<sample>::to_double(yq) - std::clamp(yd, -1.0, 1.0 - lsb)) / lsb);
                }
                std::vector<iir_state<sample>> si(tq.state_size(), iir_state<sample>{});
                std::vector<double>            sdi(td.state_size(), 0.0);
                for (std::size_t k = 0; k < x.size(); ++k) {
                    std::array<sample, 2> yq{};
                    std::array<double, 2> yd{};
                    halfband_interpolate(tq, si.data(), xq[k], yq.data());
                    halfband_interpolate(td, sdi.data(), xd[k], yd.data());
                    for (std::size_t j = 0; j < 2; ++j) {
                        track(largest_state<sample>(si),
                              std::fabs(sample_scale<sample>::to_double(yq[j]) - std::clamp(yd[j], -1.0, 1.0 - lsb))
                                  / lsb);
                    }
                }
            }
            for (const auto& d : k_near_linear_halfband_tables) {
                const near_linear_halfband_table<sample> tq(d);
                const near_linear_halfband_table<double> td = as_double(tq);
                std::vector<iir_state<sample>>           sq(tq.state_size(), iir_state<sample>{});
                std::vector<double>                      sd(td.state_size(), 0.0);
                for (std::size_t k = 0; k < x.size(); ++k) {
                    std::array<sample, 2> yq{};
                    std::array<double, 2> yd{};
                    const sample          dq = k >= tq.delay ? xq[k - tq.delay] : sample{};
                    const double          dd = k >= td.delay ? xd[k - td.delay] : 0.0;
                    near_linear_interpolate(tq, sq.data(), dq, xq[k], yq.data());
                    near_linear_interpolate(td, sd.data(), dd, xd[k], yd.data());
                    track(largest_state<sample>(sq),
                          std::fabs(sample_scale<sample>::to_double(yq[1]) - std::clamp(yd[1], -1.0, 1.0 - lsb)) / lsb);
                }
            }
            for (const auto& d : k_fourpath_tables) {
                const fourpath_table<sample>   tq(d);
                const fourpath_table<double>   td = as_double(tq);
                std::vector<iir_state<sample>> sq(tq.state_size(), iir_state<sample>{});
                std::vector<double>            sd(td.state_size(), 0.0);
                for (std::size_t k = 0; k < x.size(); ++k) {
                    std::array<sample, 4> yq{};
                    std::array<double, 4> yd{};
                    fourpath_interpolate(tq, sq.data(), xq[k], yq.data());
                    fourpath_interpolate(td, sd.data(), xd[k], yd.data());
                    for (std::size_t j = 0; j < 4; ++j) {
                        track(largest_state<sample>(sq),
                              std::fabs(sample_scale<sample>::to_double(yq[j]) - std::clamp(yd[j], -1.0, 1.0 - lsb))
                                  / lsb);
                    }
                }
            }
            for (const auto& d : k_guards) {
                const biquad_cascade_table<sample> tq(d);
                const biquad_cascade_table<double> td = as_double(tq);
                std::vector<iir_state<sample>>     sq(tq.state_size(), iir_state<sample>{});
                std::vector<double>                sd(td.state_size(), 0.0);
                for (std::size_t k = 0; k < x.size(); ++k) {
                    const sample yq = biquad_cascade(tq, sq.data(), xq[k]);
                    const double yd = biquad_cascade(td, sd.data(), xd[k]);
                    track(largest_state<sample>(sq),
                          std::fabs(sample_scale<sample>::to_double(yq) - std::clamp(yd, -1.0, 1.0 - lsb)) / lsb);
                }
            }
        }
        std::printf("[ measured ] %s: largest state %.3f of full scale (two guard bits: < 4), output within %.3f "
                    "LSB of double\n",
                    std::is_same_v<sample, std::int16_t> ? "Q15" : "Q31", worst_state, worst_err);
        EXPECT_LT(worst_state, 4.0) << "a state left the two integer guard bits";
        EXPECT_LE(worst_state, 3.0) << "the measured headroom (IIR-0: 2.30 at worst in Q15)";
        // Every output, including the clamped overshoots, within one LSB of
        // the double output (same coefficients) saturated to the format:
        // nothing wrapped inside (measured 0.5 LSB, the output's own
        // rounding, in both formats).
        EXPECT_LE(worst_err, 1.0);
    }

    TYPED_TEST(iir_fixed, OvershootSaturatesAtTheOutput) {
        // A full-scale square wave through the interpolator rings past full
        // scale (Gibbs): the output clamps there, it never wraps.
        using sample = TypeParam;
        std::vector<double> x(512);
        for (std::size_t i = 0; i < x.size(); ++i) {
            x[i] = (i / 16) % 2 == 0 ? 1.0 : -1.0;
        }
        const auto& d       = k_nonlinear[4];
        const auto  yq      = run_halfband_interpolate(nonlinear_halfband_table<sample>(d), quantize<sample>(x));
        const auto  yd      = run_halfband_interpolate(nonlinear_halfband_table<double>(d), unit(quantize<sample>(x)));
        int         clamped = 0;
        for (std::size_t i = 0; i < yq.size(); ++i) {
            if (yd[i] >= 1.0) {
                EXPECT_EQ(yq[i], std::numeric_limits<sample>::max());
                ++clamped;
            }
            else if (yd[i] < -1.0) {
                EXPECT_EQ(yq[i], std::numeric_limits<sample>::min());
                ++clamped;
            }
        }
        EXPECT_GT(clamped, 0) << "the overshoot happened";
    }

    // ------------------------------------------------------------------------
    // The limit-cycle decay (I8; I6's fixed-point window)
    // ------------------------------------------------------------------------

    /// The longest zero-input decay of one floor-truncated one-multiplier
    /// section from the full state headroom to exactly zero: the kernel's own
    /// arithmetic (x = x1 = 0), from both ends of the range (the map is
    /// monotone in |u|, so they bound every start).
    template <typename S>
    std::size_t section_decay(typename iir_traits<S>::coeff1 a) {
        using st                 = iir_state<S>;
        std::size_t        worst = 0;
        const std::int64_t top   = std::int64_t{1} << (iir_traits<S>::k_state_frac_bits + 2);
        for (const std::int64_t start : {top - 1, -top}) {
            auto        u = static_cast<st>(start);
            std::size_t n = 0;
            while (u != 0 && n < 100000) {
                u = iir_traits<S>::allpass1(st{}, st{}, u, a);
                ++n;
            }
            worst = std::max(worst, n);
        }
        return worst;
    }

    TYPED_TEST(iir_fixed, SilenceDecaysToAnExactlyZeroStateWithinTheBound) {
        using sample                  = TypeParam;
        const auto  inputs            = full_scale_inputs(1200);
        std::size_t worst_nonlinear   = 0;
        std::size_t worst_bound       = 0;
        std::size_t worst_near_linear = 0;
        std::size_t worst_guard       = 0;
        for (const auto& d : k_nonlinear) {
            const nonlinear_halfband_table<sample> t(d);
            // IIR-0's bound: per path, 1 + the sum of its sections' decays;
            // the stage, the longer path plus the decimator's odd latch.
            std::array<std::size_t, 2> path{1, 1};
            for (std::size_t i = 0; i < t.size0; ++i) {
                path[0] += section_decay<sample>(t.path0[i]);
            }
            for (std::size_t i = 0; i < t.size1; ++i) {
                path[1] += section_decay<sample>(t.path1[i]);
            }
            const std::size_t bound = std::max(path[0], path[1]) + 1;
            worst_bound             = std::max(worst_bound, bound);
            for (const auto& x : inputs) {
                const std::vector<sample>      xq = quantize<sample>(x);
                std::vector<iir_state<sample>> st(t.state_size(), iir_state<sample>{});
                sample                         odd{};
                for (std::size_t k = 0; k + 1 < xq.size(); k += 2) {
                    static_cast<void>(halfband_decimate(t, st.data(), xq[k], odd));
                    odd = xq[k + 1];
                }
                std::size_t frames = 0;
                while (std::any_of(st.begin(), st.end(), [](auto v) { return v != 0; }) && frames < 10 * bound) {
                    static_cast<void>(halfband_decimate(t, st.data(), sample{}, odd));
                    odd = sample{};
                    ++frames;
                }
                EXPECT_LE(frames, bound) << "order " << d.order;
                EXPECT_EQ(halfband_decimate(t, st.data(), sample{}, sample{}), sample{}) << "and the output is silence";
                worst_nonlinear = std::max(worst_nonlinear, frames);
            }
        }
        for (const auto& d : k_near_linear_halfband_tables) {
            const near_linear_halfband_table<sample> t(d);
            for (const auto& x : inputs) {
                const std::vector<sample>      xq = quantize<sample>(x);
                std::vector<iir_state<sample>> st(t.state_size(), iir_state<sample>{});
                for (std::size_t k = 0; k < xq.size(); ++k) {
                    std::array<sample, 2> y{};
                    near_linear_interpolate(t, st.data(), sample{}, xq[k], y.data());
                }
                std::size_t frames = 0;
                while (std::any_of(st.begin(), st.end(), [](auto v) { return v != 0; }) && frames < 100000) {
                    std::array<sample, 2> y{};
                    near_linear_interpolate(t, st.data(), sample{}, sample{}, y.data());
                    ++frames;
                }
                worst_near_linear = std::max(worst_near_linear, frames);
            }
        }
        for (const auto& d : k_guards) {
            const biquad_cascade_table<sample> t(d);
            for (const auto& x : inputs) {
                const std::vector<sample>      xq = quantize<sample>(x);
                std::vector<iir_state<sample>> st(t.state_size(), iir_state<sample>{});
                for (const sample v : xq) {
                    static_cast<void>(biquad_cascade(t, st.data(), v));
                }
                std::size_t frames = 0;
                while (std::any_of(st.begin(), st.end(), [](auto v) { return v != 0; }) && frames < 100000) {
                    static_cast<void>(biquad_cascade(t, st.data(), sample{}));
                    ++frames;
                }
                worst_guard = std::max(worst_guard, frames);
            }
        }
        std::printf("[ measured ] %s: silence to an exactly zero state within %zu frames (nonlinear; bound %zu), %zu "
                    "(near-linear), %zu (guard band cascade)\n",
                    std::is_same_v<sample, std::int16_t> ? "Q15" : "Q31", worst_nonlinear, worst_bound,
                    worst_near_linear, worst_guard);
        // IIR-0's magnitude-truncation bounds at the largest designs (Q15 / Q31,
        // 120 dB, 5/12, d = 1): near-linear 518 / 1003 low-rate frames.
        EXPECT_LE(worst_near_linear, (std::is_same_v<sample, std::int16_t> ? 518u : 1003u));
        EXPECT_LT(worst_guard, 100000u) << "no zero-input limit cycle in the guard band cascade";
    }

    // ------------------------------------------------------------------------
    // The Q15 floor (IIR-0's measurement)
    // ------------------------------------------------------------------------

    TEST(IirFixedFloor, Q15FloorIsTheOutputsOwnRounding) {
        // A sine at 997/48000 of the higher rate, -1 dBFS, quantized to Q15,
        // through the Q15 decimator and through the double structure on the
        // same input; error power re a full-scale sine's, 2^14 outputs after
        // 4096. Beside it the double output merely rounded to Q15.
        for (const std::size_t i : {0u, 4u, 8u, 12u}) {
            const auto&         d = k_nonlinear[i];
            const std::size_t   n = 2 * (4096 + 16384);
            std::vector<double> x(n);
            for (std::size_t k = 0; k < n; ++k) {
                x[k] = std::pow(10.0, -1.0 / 20.0)
                       * std::sin(2.0 * std::numbers::pi * 997.0 / 48000.0 * static_cast<double>(k));
            }
            const auto xq = quantize<std::int16_t>(x);
            const auto yq = run_halfband_decimate(nonlinear_halfband_table<std::int16_t>(d), xq);
            const auto yd = run_halfband_decimate(nonlinear_halfband_table<double>(d), unit(xq));
            double     eq = 0.0;
            double     er = 0.0;
            for (std::size_t k = 4096; k < yq.size(); ++k) {
                const double q = sample_scale<std::int16_t>::to_double(yq[k]);
                const double r = sample_scale<std::int16_t>::to_double(sample_scale<std::int16_t>::from_double(yd[k]));
                eq += (q - yd[k]) * (q - yd[k]);
                er += (r - yd[k]) * (r - yd[k]);
            }
            const double m   = static_cast<double>(yq.size() - 4096);
            const double edb = 10.0 * std::log10(eq / m / 0.5);
            const double rdb = 10.0 * std::log10(er / m / 0.5);
            std::printf("[ measured ] Q15 floor, order %d: %.2f dB (the double output rounded to Q15: %.2f dB)\n",
                        d.order, edb, rdb);
            EXPECT_NEAR(edb, rdb, 0.1);
            EXPECT_LE(edb, -97.5);
        }
    }

    // ------------------------------------------------------------------------
    // The bits
    // ------------------------------------------------------------------------

    template <typename S>
    std::vector<std::uint64_t> fingerprints() {
        const std::vector<double> xd = [] {
            xorshift32          rng(0x5eedu);
            std::vector<double> v(1024);
            for (double& s : v) {
                s = 0.9 * rng.next_unit();
            }
            return v;
        }();
        const std::vector<S>       x = quantize<S>(xd);
        std::vector<std::uint64_t> out;
        std::vector<S>             all;
        const auto                 append = [&](const std::vector<S>& y) { all.insert(all.end(), y.begin(), y.end()); };
        for (const auto& d : k_nonlinear) {
            append(run_halfband_decimate(nonlinear_halfband_table<S>(d), x));
        }
        out.push_back(fnv1a64(all));
        all.clear();
        for (const auto& d : k_nonlinear) {
            append(run_halfband_interpolate(nonlinear_halfband_table<S>(d), x));
        }
        out.push_back(fnv1a64(all));
        all.clear();
        for (const auto& d : k_near_linear_halfband_tables) {
            append(run_near_linear_decimate(near_linear_halfband_table<S>(d), x));
            append(run_near_linear_interpolate(near_linear_halfband_table<S>(d), x));
        }
        out.push_back(fnv1a64(all));
        all.clear();
        for (const auto& d : k_fourpath_tables) {
            append(run_fourpath_decimate(fourpath_table<S>(d), x));
            append(run_fourpath_interpolate(fourpath_table<S>(d), x));
        }
        out.push_back(fnv1a64(all));
        all.clear();
        for (const auto& d : k_guards) {
            append(run_biquad_cascade(biquad_cascade_table<S>(d), x));
        }
        out.push_back(fnv1a64(all));
        return out;
    }

    TYPED_TEST(iir_fixed, OutputBitsMatchThePins) {
        using sample    = TypeParam;
        const auto  got = fingerprints<sample>();
        const auto& pin = std::is_same_v<sample, std::int16_t> ? tap::dsp::test::iir_fixed_pins::k_q15
                                                               : tap::dsp::test::iir_fixed_pins::k_q31;
        ASSERT_EQ(got.size(), std::size(pin));
        for (std::size_t i = 0; i < got.size(); ++i) {
            std::printf("[ measured ] %s %-30s 0x%016" PRIx64 "\n",
                        std::is_same_v<sample, std::int16_t> ? "Q15" : "Q31", pin[i].name, got[i]);
            EXPECT_EQ(got[i], pin[i].fingerprint) << pin[i].name << ": a fixed-point output bit moved";
        }
    }

} // namespace
