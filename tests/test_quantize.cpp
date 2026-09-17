// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Contract battery for row-sum-preserving quantization, ported from the
// row-sum checks in SampleRateTap's test_fixed_point.cpp (there stated
// against the assembled polyphase bank; here against the utility itself,
// over rows of a genuinely designed prototype).

#include <cmath>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/kaiser.h"
#include "tap/dsp/quantize.h"

namespace {

    using tap::dsp::quantize_row_preserving_sum;
    using tap::dsp::sample_traits;

    // A realistic source: the RatioTap-shaped design (L = 160, plain Kaiser)
    // whose branches this utility will quantize in production.
    std::vector<double> designed_prototype(std::size_t phases, std::size_t taps) {
        std::vector<double> h(phases * taps);
        tap::dsp::design_prototype(h, phases, (19000.0 + 22050.0) / 48000.0, tap::dsp::kaiser_beta(70.0));
        return h;
    }

    template <typename S>
    void check_rows_sum_exact() {
        constexpr std::size_t k_phases = 160;
        constexpr std::size_t k_taps   = 32;
        const auto            proto    = designed_prototype(k_phases, k_taps);
        using coeff                    = typename sample_traits<S>::coeff;
        const double scale             = sample_traits<S>::k_coeff_scale;

        std::vector<double> row(k_taps);
        std::vector<coeff>  q(k_taps);
        for (std::size_t p = 0; p < k_phases; ++p) {
            double exact = 0.0;
            for (std::size_t t = 0; t < k_taps; ++t) {
                row[t] = proto[t * k_phases + p];
                exact += row[t] * scale;
            }
            quantize_row_preserving_sum<S>(row, q);
            std::int64_t sum = 0;
            for (std::size_t t = 0; t < k_taps; ++t) {
                sum += q[t];
                // Each tap stays within the rounding step plus at most the
                // residual corrections: never grossly redistributed.
                EXPECT_LE(std::abs(static_cast<double>(q[t]) - row[t] * scale), 1.5) << "phase " << p << " tap " << t;
            }
            ASSERT_EQ(sum, std::llround(exact)) << "phase " << p;
        }
    }

    TEST(Quantize, RowSumsAreExactQ15) {
        check_rows_sum_exact<std::int16_t>();
    }

    TEST(Quantize, RowSumsAreExactQ31) {
        check_rows_sum_exact<std::int32_t>();
    }

    TEST(Quantize, DoubleIsPlainConversion) {
        // The golden model's coefficients are the designed values themselves.
        const std::vector<double> row{0.25, -0.125, 1.0, -0.9999, 0.0, 1e-300};
        std::vector<double>       q(row.size());
        quantize_row_preserving_sum<double>(row, q);
        EXPECT_EQ(q, row);
    }

    // The +/-1 correction never wraps a tap at the rail. A tap that saturated
    // in make_coeff has the row's largest positive remainder, so a naive
    // largest-remainder pick would step it: instead the step goes to the
    // largest remainder among taps that can still move, and the sum is
    // preserved through them. When every tap sits at the needed rail the
    // correction stops and the row keeps its saturated sum.
    TEST(Quantize, CorrectionNeverWrapsATapAtTheRail) {
        // Just over the rail: 2.0001 * 16384 = 32769.6 -> 32767 (3 LSB short).
        const std::vector<double> row{2.0001, 0.4, -0.25};
        std::vector<std::int16_t> q15(row.size());
        quantize_row_preserving_sum<std::int16_t>(row, q15);
        EXPECT_EQ(q15[0], 32767); // saturated, not wrapped, untouched
        EXPECT_EQ(std::int64_t{q15[0]} + q15[1] + q15[2], std::llround((2.0001 + 0.4 - 0.25) * 16384));
        std::vector<std::int32_t> q31(row.size());
        quantize_row_preserving_sum<std::int32_t>(row, q31);
        EXPECT_EQ(q31[0], 2147483647);
        EXPECT_EQ(std::int64_t{q31[0]} + q31[1] + q31[2], std::llround((2.0001 + 0.4 - 0.25) * 1073741824.0));
        // The negative rail too.
        const std::vector<double> neg{-2.0001, -0.4};
        std::vector<std::int16_t> n15(neg.size());
        quantize_row_preserving_sum<std::int16_t>(neg, n15);
        EXPECT_EQ(n15[0], -32768);
        EXPECT_EQ(std::int64_t{n15[0]} + n15[1], std::llround((-2.0001 - 0.4) * 16384));
        // Every tap at the rail: nothing can move, the loop stops, no wrap.
        const std::vector<double> all{2.5, 3.0};
        std::vector<std::int16_t> a15(all.size());
        quantize_row_preserving_sum<std::int16_t>(all, a15);
        EXPECT_EQ(a15[0], 32767);
        EXPECT_EQ(a15[1], 32767);
        const std::vector<double> alln{-2.5, -3.0};
        quantize_row_preserving_sum<std::int16_t>(alln, a15);
        EXPECT_EQ(a15[0], -32768);
        EXPECT_EQ(a15[1], -32768);
    }

    // The rule that makes the rail case honest: a saturated tap has the
    // row's largest positive remainder, so a naive largest-remainder pick
    // would land on it and either wrap (main before this) or give up on a
    // sum that the other taps can still absorb. 1.99997 is the first value
    // that saturates in Q1.14 (32767.5 rounds to the rail); the row's sum
    // 42597.9 -> 42598 is representable and must be reached through taps 1
    // and 2.
    TEST(Quantize, NearRailRowStillPreservesItsSum) {
        const std::vector<double> row{1.99997, 0.3, 0.3};
        std::vector<std::int16_t> q15(row.size());
        quantize_row_preserving_sum<std::int16_t>(row, q15);
        EXPECT_EQ(q15[0], 32767);                    // at the rail, untouched
        EXPECT_EQ(q15[1] + q15[2], 4915 + 4915 + 1); // the step went to a movable tap
        EXPECT_EQ(std::int64_t{q15[0]} + q15[1] + q15[2], std::llround(1.99997 * 16384 + 0.6 * 16384));
        const std::vector<double> row31{1.9999999995, 0.3, 0.3};
        std::vector<std::int32_t> q31(row31.size());
        quantize_row_preserving_sum<std::int32_t>(row31, q31);
        EXPECT_EQ(q31[0], 2147483647);
        EXPECT_EQ(std::int64_t{q31[0]} + q31[1] + q31[2], std::llround((1.9999999995 + 0.6) * 1073741824.0));
        EXPECT_EQ(std::int64_t{q31[0]} + q31[1] + q31[2], 2791728742);
    }

    TEST(Quantize, FloatIsPlainConversion) {
        const std::vector<double> row{0.25, -0.125, 1.0, -0.9999, 0.0};
        std::vector<float>        q(row.size());
        quantize_row_preserving_sum<float>(row, q);
        for (std::size_t t = 0; t < row.size(); ++t) {
            EXPECT_FLOAT_EQ(q[t], static_cast<float>(row[t]));
        }
    }

    // The consequence the correction exists for: a unity-DC row quantizes to a
    // row whose coefficient sum is exactly the format's unity, so DC gain
    // through any phase deviates by at most one output LSB.
    TEST(Quantize, UnityDcRowSumsToFormatUnity) {
        constexpr std::size_t k_taps = 48;
        // A smooth row normalized to sum exactly 1.0 in double.
        std::vector<double> row(k_taps);
        double              sum = 0.0;
        for (std::size_t t = 0; t < k_taps; ++t) {
            const double u = (static_cast<double>(t) - 23.5) / 24.0;
            row[t]         = std::exp(-4.0 * u * u);
            sum += row[t];
        }
        for (auto& v : row) {
            v /= sum;
        }
        std::vector<std::int16_t> q(k_taps);
        quantize_row_preserving_sum<std::int16_t>(row, q);
        std::int64_t qsum = 0;
        for (const auto c : q) {
            qsum += c;
        }
        EXPECT_EQ(qsum, 16384); // Q1.14 unity, exactly
    }

} // namespace
