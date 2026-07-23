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
