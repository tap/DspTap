// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// butterworth.h without the host's flush-to-zero: the honest limit
// (subnormals.h, IIR_PLAN I9). Built into tap_dsp_iir_subnormals with
// TAP_DSP_ALLOW_SUBNORMALS, which compiles the kernel's debug-build FTZ
// precondition out (tests/CMakeLists.txt says why it is an executable of its
// own). The tail after window_frames() stays below the same bound as under
// FTZ (test_butterworth.cpp), and below the smallest normal value: without
// FTZ the residue is subnormal, never the limit cycle FTZ can hold.

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "support/iir_support.h"
#include "tap/dsp/butterworth.h"
#include "tap/dsp/subnormals.h"

namespace {

    using tap::dsp::basic_butterworth_lowpass;
    using tap::dsp::test::no_flush_to_zero_fixture;

    class butterworth_subnormals : public no_flush_to_zero_fixture {};

    template <typename S>
    std::vector<S> square(std::size_t n, std::size_t half_period, S amp = S(1)) {
        std::vector<S> x(n);
        for (std::size_t i = 0; i < n; ++i) {
            x[i] = ((i / half_period) % 2) != 0 ? amp : -amp;
        }
        return x;
    }

    /// The largest |output| after the window, from a full-scale square drive
    /// and from one near the format's largest value.
    template <typename S, std::size_t N>
    S tail_after_window(double fc) {
        basic_butterworth_lowpass<S, N> f(fc);
        S                               worst = S(0);
        for (const S amp : {S(1), static_cast<S>(std::numeric_limits<S>::max() / S(1e6))}) {
            f.reset();
            auto x = square<S>(4000, 5, amp);
            f.process(x.data(), x.size(), x.data());
            std::vector<S> z(2 * f.window_frames(), S(0));
            f.process(z.data(), z.size(), z.data());
            for (std::size_t i = f.window_frames(); i < z.size(); ++i) {
                worst = std::max(worst, static_cast<S>(std::fabs(z[i])));
            }
        }
        return worst;
    }

    TEST_F(butterworth_subnormals, TheTailAfterTheWindowIsBounded) {
        // float at every cutoff the FTZ battery runs; double (software on the
        // M4, M4F and M33 legs) at two.
        for (const double fc : {0.005, 0.02, 0.1, 0.25, 0.45}) {
            EXPECT_LT((tail_after_window<float, 4>(fc)), std::numeric_limits<float>::min()) << fc;
            EXPECT_LT((tail_after_window<float, 8>(fc)), std::numeric_limits<float>::min()) << fc;
        }
        for (const double fc : {0.1, 0.25}) {
            EXPECT_LT((tail_after_window<double, 4>(fc)), std::numeric_limits<double>::min()) << fc;
        }
    }

    // The design that holds a limit cycle under FTZ (test_butterworth.cpp's
    // UnderFlushToZeroSomeDesignsNeverFallSilent) decays below the smallest
    // normal value without it.
    TEST_F(butterworth_subnormals, WithoutFlushToZeroTheResidueIsSubnormal) {
        basic_butterworth_lowpass<float, 8> f(0.45);
        auto                                x = square<float>(4000, 5);
        f.process(x.data(), x.size(), x.data());
        std::vector<float> z(100000, 0.0f);
        f.process(z.data(), z.size(), z.data());
        float tail = 0.0f;
        for (std::size_t i = z.size() - 1000; i < z.size(); ++i) {
            tail = std::max(tail, std::fabs(z[i]));
        }
        EXPECT_LT(tail, std::numeric_limits<float>::min());
    }

    TEST_F(butterworth_subnormals, NeverWritesTheFloatingPointEnvironment) {
        using tap::dsp::test::k_control_mask;
        const auto                          before = tap::dsp::floating_point_control() & k_control_mask;
        basic_butterworth_lowpass<float, 4> a(0.1);
        auto                                x = square<float>(512, 3);
        a.process(x.data(), x.size(), x.data());
        a.reset();
        EXPECT_EQ(tap::dsp::floating_point_control() & k_control_mask, before);
    }

} // namespace
