// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Contract battery for row-sum-preserving quantization, ported from the
// row-sum checks in SampleRateTap's test_fixed_point.cpp (there stated
// against the assembled polyphase bank; here against the utility itself,
// over rows of a genuinely designed prototype). The symmetric form
// (quantize_symmetric_row_preserving_sum) is pinned over every committed
// guard band row (iir_tables.h: the double rows and the Q15 rows) and a
// Kaiser custom design, in all four formats: its symmetry, its zero at
// Nyquist, its sum, and its bit-identity with the half-then-mirror procedure
// SampleRateTap's guard band filter ran inline.

#include <cmath>
#include <cstdint>
#include <cstring>
#include <span>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/iir_design.h"
#include "tap/dsp/kaiser.h"
#include "tap/dsp/quantize.h"

namespace {

    using tap::dsp::quantize_row_preserving_sum;
    using tap::dsp::quantize_symmetric_row_preserving_sum;
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

    // The tie rule: equal remainders step the lowest index. A mirrored row
    // (every designer's output) therefore quantizes identically on every
    // host, since its mirrored pairs tie exactly; and the one step a tie
    // absorbs cannot be split, so the quantized row is NOT mirrored.
    TEST(Quantize, TiesStepTheLowestIndex) {
        // Scaled: 4915.45, 6553.2, 4915.45 -> rounded 4915, 6553, 4915 (sum
        // 16383) against the exact 16384.1 -> 16384: one step up, and the
        // outer taps tie for the largest remainder (0.45 against 0.2).
        const std::vector<double> row{0.3 + 0.25 / 16384.0, 0.4 - 0.4 / 16384.0, 0.3 + 0.25 / 16384.0};
        std::vector<std::int16_t> q(row.size());
        quantize_row_preserving_sum<std::int16_t>(row, q);
        EXPECT_EQ(std::int64_t{q[0]} + q[1] + q[2], 16384);
        EXPECT_EQ(q[0], 4916); // the tie's lower index took the step
        EXPECT_EQ(q[1], 6553);
        EXPECT_EQ(q[2], 4915);
        // The row's mirror is the row, and it quantizes to the same table:
        // the step lands at index 0 again, not at the mirrored index.
        const std::vector<double> mirror{row[2], row[1], row[0]};
        std::vector<std::int16_t> qm(mirror.size());
        quantize_row_preserving_sum<std::int16_t>(mirror, qm);
        EXPECT_EQ(qm, q);
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

    // ------------------------------------------------------------------------
    // The symmetric (type II) form
    // ------------------------------------------------------------------------

    /// Every even-length symmetric row the family quantizes this way: the
    /// committed guard band FIRs (double and Q15 rows) and a Kaiser custom
    /// design at 70 dB, 3/8, the order-5 Chebyshev II's edge (52 taps), as
    /// SampleRateTap's guard band filter designs a custom setting.
    std::vector<std::vector<double>> symmetric_rows() {
        std::vector<std::vector<double>> rows;
        for (const auto& r : tap::dsp::k_guard_band_fir_tables) {
            rows.emplace_back(r.taps.begin(), r.taps.end());
        }
        for (const auto& r : tap::dsp::k_guard_band_fir_q15_tables) {
            rows.emplace_back(r.taps.begin(), r.taps.end());
        }
        const double p    = 3.0 / 8.0;
        const double edge = tap::dsp::design_cheby2_lowpass(5, p, 70).stop_edge;
        rows.emplace_back(52);
        tap::dsp::design_guard_band_linear(rows.back(), p, edge, 70);
        return rows;
    }

    /// SampleRateTap's guard band filter's build_rows() (rational's
    /// guard_band.h, IIR-c), as it ran before this helper, its members made
    /// locals: the first half quantized into a scratch row of N / 2, copied
    /// into both halves.
    template <typename S>
    std::vector<typename sample_traits<S>::coeff> half_then_mirror(const std::vector<double>& design) {
        using coeff             = typename sample_traits<S>::coeff;
        const std::size_t  n    = design.size();
        const std::size_t  half = n / 2;
        std::vector<coeff> q(half);
        quantize_row_preserving_sum<S>(std::span<const double>(design.data(), half), q);
        std::vector<coeff> rows;
        rows.assign(n, coeff{});
        for (std::size_t i = 0; i < half; ++i) {
            rows[i]         = q[i];
            rows[n - 1 - i] = q[i];
        }
        return rows;
    }

    template <typename S>
    std::vector<typename sample_traits<S>::coeff> symmetric(const std::vector<double>& design) {
        std::vector<typename sample_traits<S>::coeff> q(design.size());
        quantize_symmetric_row_preserving_sum<S>(design, q);
        return q;
    }

    template <typename S>
    void check_half_then_mirror_identity() {
        for (const auto& design : symmetric_rows()) {
            const auto a = symmetric<S>(design);
            const auto b = half_then_mirror<S>(design);
            ASSERT_EQ(a.size(), b.size());
            EXPECT_EQ(std::memcmp(a.data(), b.data(), a.size() * sizeof(a[0])), 0)
                << design.size() << " taps: bit for bit";
        }
    }

    TEST(Quantize, SymmetricRowIsTheHalfQuantizedThenMirroredBitForBit) {
        // The identity SampleRateTap's switch to this helper rests on: its
        // guard band filter's table, in every format, is unchanged by it.
        check_half_then_mirror_identity<double>();
        check_half_then_mirror_identity<float>();
        check_half_then_mirror_identity<std::int16_t>();
        check_half_then_mirror_identity<std::int32_t>();
    }

    template <typename S>
    void check_symmetric_with_zero_at_nyquist() {
        for (const auto& design : symmetric_rows()) {
            const auto        q = symmetric<S>(design);
            const std::size_t n = q.size();
            for (std::size_t i = 0; i < n / 2; ++i) {
                EXPECT_EQ(std::memcmp(&q[i], &q[n - 1 - i], sizeof(q[i])), 0) << "symmetric, bit for bit";
            }
            if constexpr (sample_traits<S>::k_is_fixed_point) {
                // H(1/2) = sum_i (-1)^i q[i], in integers: exactly zero.
                std::int64_t alt = 0;
                for (std::size_t i = 0; i < n; ++i) {
                    alt += (i % 2 == 0 ? 1 : -1) * static_cast<std::int64_t>(q[i]);
                }
                EXPECT_EQ(alt, 0) << n << " taps: the zero at Nyquist, exactly, in the quantized taps";
            }
            else {
                // Each mirrored pair cancels exactly.
                double alt = 0.0;
                for (std::size_t i = 0; i < n / 2; ++i) {
                    alt += (i % 2 == 0 ? 1.0 : -1.0) * (static_cast<double>(q[i]) - static_cast<double>(q[n - 1 - i]));
                }
                EXPECT_EQ(alt, 0.0);
            }
        }
    }

    TEST(Quantize, SymmetricRowIsSymmetricWithAnExactZeroAtNyquist) {
        check_symmetric_with_zero_at_nyquist<double>();
        check_symmetric_with_zero_at_nyquist<float>();
        check_symmetric_with_zero_at_nyquist<std::int16_t>();
        check_symmetric_with_zero_at_nyquist<std::int32_t>();
    }

    template <typename S>
    void check_symmetric_sum() {
        const double scale = sample_traits<S>::k_coeff_scale;
        for (const auto& design : symmetric_rows()) {
            const auto        q    = symmetric<S>(design);
            const std::size_t n    = q.size();
            double            half = 0.0; // the half's scaled sum, as quantize_row_preserving_sum forms it
            double            full = 0.0; // the design's
            for (std::size_t t = 0; t < n; ++t) {
                if (t < n / 2) {
                    half += design[t] * scale;
                }
                full += design[t] * scale;
            }
            std::int64_t sum = 0;
            for (const auto c : q) {
                sum += c;
            }
            EXPECT_EQ(sum, 2 * std::llround(half)) << n << " taps: twice the half's preserved sum";
            EXPECT_LE(std::fabs(static_cast<double>(sum) - full), 1.0) << "within one LSB of the design's sum";
        }
    }

    TEST(Quantize, SymmetricRowSumIsTwiceTheHalfsPreservedSum) {
        check_symmetric_sum<std::int16_t>();
        check_symmetric_sum<std::int32_t>();
    }

    TEST(Quantize, SymmetricRowStaysSymmetricWhereTheWholeRowWouldNot) {
        // {a, b, b, a} scaled to Q1.14: 4915.45, 3276.85, 3276.85, 4915.45,
        // rounding to 4915 and 3277. The whole row's sum 16384.6 wants 16385,
        // one step, which the tied outer pair cannot split: the step goes to
        // index 0, and the row loses its symmetry and its zero at Nyquist.
        // The half's sum 8192.3 wants 8192, which it has: no step, and the
        // mirrored row sums to 16384, within one LSB of 16384.6.
        const double              a = 4915.45 / 16384.0;
        const double              b = 3276.85 / 16384.0;
        const std::vector<double> row{a, b, b, a};
        std::vector<std::int16_t> whole(row.size());
        quantize_row_preserving_sum<std::int16_t>(row, whole);
        EXPECT_EQ(whole, (std::vector<std::int16_t>{4916, 3277, 3277, 4915}));
        EXPECT_NE(whole[0] - whole[1] + whole[2] - whole[3], 0) << "the whole row's Nyquist zero is lost";
        std::vector<std::int16_t> sym(row.size());
        quantize_symmetric_row_preserving_sum<std::int16_t>(row, sym);
        EXPECT_EQ(sym, (std::vector<std::int16_t>{4915, 3277, 3277, 4915}));
        EXPECT_EQ(sym[0] - sym[1] + sym[2] - sym[3], 0);
    }

    TEST(Quantize, SymmetricRowIsPlainConversionInTheFloatingFormats) {
        // No correction runs (k_coeff_scale 1): the first half converted,
        // then mirrored, which for a symmetric design is the whole row
        // converted.
        for (const auto& design : symmetric_rows()) {
            const auto d = symmetric<double>(design);
            const auto f = symmetric<float>(design);
            for (std::size_t t = 0; t < design.size(); ++t) {
                EXPECT_EQ(d[t], design[t]);
                EXPECT_EQ(f[t], static_cast<float>(design[t]));
            }
        }
    }

} // namespace
