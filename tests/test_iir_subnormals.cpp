// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The honest limit of the recursive kernels without the host's flush-to-zero
// (tap/dsp/subnormals.h, IIR_PLAN I9), in an executable of its own built with
// TAP_DSP_ALLOW_SUBNORMALS (tests/CMakeLists.txt: the macro is an ODR rule,
// one setting per image, so it cannot share tap_dsp_tests):
//   - WithoutFlushToZeroTheResidueIsBounded: with the bits cleared (or on a
//     target where none governs the type: double on the M4F and M33,
//     everything on the soft-float M4), silence after signal leaves the
//     states at a residue of a few multiples of the smallest subnormal, which
//     the stated bound holds, and the output there is that residue, not
//     silence; the counterpart of iir_kernel.SilenceDecaysToExactZero;
//   - NeverTouchesTheFloatingPointEnvironment in the other host state (the
//     bits clear), the counterpart of IirKernel's;
//   - the precondition compiled out: in a Debug build these calls would
//     assert in tap_dsp_tests (IirKernelDebug), and here they run.
// Every case runs on the four QEMU legs.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "support/iir_support.h"
#include "support/signals.h"
#include "tap/dsp/iir_design.h"
#include "tap/dsp/iir_kernels.h"
#include "tap/dsp/subnormals.h"

#if !defined(TAP_DSP_ALLOW_SUBNORMALS)
#error "tap_dsp_iir_subnormals is built with TAP_DSP_ALLOW_SUBNORMALS (tests/CMakeLists.txt)"
#endif

namespace {

    using namespace tap::dsp;
    using namespace tap::dsp::test;

    template <typename T>
    class iir_subnormals : public no_flush_to_zero_fixture {};

    using floating_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(iir_subnormals, floating_types);

    /// The honest limit's bound, in multiples of the smallest subnormal: a
    /// first-order section with coefficient a sticks at k of them for k <=
    /// 1 / (2 (1 - a)) (under round to nearest a k rounds back to k); a path
    /// of n sections, each holding at most that, and the output, a half-sum
    /// of two such paths, at most n times it. Stated as one number per design.
    std::size_t residue_bound(const nonlinear_halfband_design& d) {
        const double a = d.coefficients[d.count - 1];
        const auto   k = static_cast<std::size_t>(std::floor(1.0 / (2.0 * (1.0 - a))));
        return (k + 1) * (d.count + 2);
    }

    TYPED_TEST(iir_subnormals, WithoutFlushToZeroTheResidueIsBounded) {
        using real = TypeParam;
        EXPECT_FALSE(subnormals_flushed<real>()) << "the fixture cleared the bits, or none governs real here";
        constexpr double dmin = static_cast<double>(std::numeric_limits<real>::denorm_min());
        xorshift32       rng(0x7777u);
        for (const auto& d : {k_nonlinear_halfband<70, 3, 8, 1>, k_nonlinear_halfband<120, 5, 12, 1>,
                              k_nonlinear_halfband<70, 1, 3, 2>}) {
            const nonlinear_halfband_table<real> t(d);
            std::vector<iir_state<real>>         st(t.state_size(), iir_state<real>{});
            for (int k = 0; k < 2000; ++k) {
                static_cast<void>(halfband_decimate(t, st.data(), static_cast<real>(rng.next_unit()),
                                                    static_cast<real>(rng.next_unit())));
            }
            const std::size_t window = floating_flush_frames<real>(d.coefficients[d.count - 1], t.multiplies());
            double            out    = 0.0;
            for (std::size_t k = 0; k < window + 4096; ++k) {
                const real y = halfband_decimate(t, st.data(), real(0), real(0));
                if (k >= window) {
                    out = std::max(out, std::fabs(static_cast<double>(y)));
                }
            }
            double state = 0.0;
            for (const auto s : st) {
                state = std::max(state, std::fabs(static_cast<double>(s)));
            }
            const std::size_t bound = residue_bound(d);
            std::printf("[ measured ] %s, order %d, a = %.4f: after the window the state holds %.0f x and the output "
                        "%.0f x the smallest subnormal (bound %lu)\n",
                        std::is_same_v<real, float> ? "float" : "double", d.order, d.coefficients[d.count - 1],
                        state / dmin, out / dmin, static_cast<unsigned long>(bound));
            EXPECT_LE(state, static_cast<double>(bound) * dmin);
            EXPECT_LE(out, static_cast<double>(bound) * dmin) << "a residue, not a signal";
            EXPECT_LT(state, static_cast<double>(std::numeric_limits<real>::min())) << "only subnormals are left";
        }
    }

    TYPED_TEST(iir_subnormals, TheGuardBandCascadesResidueIsBounded) {
        // Without flush-to-zero the floating biquad cascade's state settles at
        // a zero-input limit cycle of a few times the smallest subnormal (the
        // round-to-nearest cycle of fixed point, at the subnormal LSB):
        // measured 6 to 20 times in float and double, every order.
        using real            = TypeParam;
        constexpr double dmin = static_cast<double>(std::numeric_limits<real>::denorm_min());
        for (const auto& d : {k_guard_band_cheby2<70, 1, 3, 5>, k_guard_band_cheby2<70, 3, 8, 7>,
                              k_guard_band_cheby2<120, 5, 12, 11>}) {
            const biquad_cascade_table<real> t(d);
            std::vector<iir_state<real>>     st(t.state_size(), iir_state<real>{});
            for (int k = 0; k < 3000; ++k) {
                static_cast<void>(biquad_cascade(t, st.data(), (k % 10) < 5 ? real(1) : real(-1)));
            }
            const std::size_t window = floating_flush_frames<real>(d.pole_radius[d.count - 1], d.count);
            double            out    = 0.0;
            for (std::size_t n = 0; n < 3 * window; ++n) {
                const real y = biquad_cascade(t, st.data(), real(0));
                if (n >= 2 * window) {
                    out = std::max(out, std::fabs(static_cast<double>(y)));
                }
            }
            double state = 0.0;
            for (const auto s : st) {
                state = std::max(state, std::fabs(static_cast<double>(s)));
            }
            std::printf("[ measured ] %s Chebyshev II order %d without flush-to-zero: the state holds %.0f x and the "
                        "output %.0f x the smallest subnormal\n",
                        std::is_same_v<real, float> ? "float" : "double", d.order, state / dmin, out / dmin);
            EXPECT_LE(state, 32.0 * dmin);
            EXPECT_LE(out, 32.0 * dmin);
        }
    }

    TYPED_TEST(iir_subnormals, NeverTouchesTheFloatingPointEnvironment) {
        using real = TypeParam;
        const nonlinear_halfband_table<real>   hb(k_nonlinear_halfband<70, 3, 8, 1>);
        const near_linear_halfband_table<real> nl(k_near_linear_halfband<70, 3, 8, 1>);
        const biquad_cascade_table<real>       bq(k_guard_band_cheby2<70, 3, 8, 7>);
        std::vector<iir_state<real>>           s1(hb.state_size()), s2(nl.state_size()), s3(bq.state_size());
        xorshift32                             rng(0x8888u);
        std::array<real, 2>                    out{};
        for (int i = 0; i < 64; ++i) {
            const real          x      = static_cast<real>(rng.next_unit());
            const std::uint64_t before = floating_point_control() & k_control_mask;
            static_cast<void>(halfband_decimate(hb, s1.data(), x, x));
            halfband_interpolate(hb, s1.data(), x, out.data());
            static_cast<void>(near_linear_decimate(nl, s2.data(), x, x));
            static_cast<void>(biquad_cascade(bq, s3.data(), x));
            EXPECT_EQ(floating_point_control() & k_control_mask, before);
        }
    }

} // namespace
