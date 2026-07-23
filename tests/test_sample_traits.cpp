// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Contract battery for the sample-format traits, ported and extended from
// SampleRateTap's test_fixed_point.cpp. Every number here is a documented
// contract point (Q formats, rounding mode, saturation, accumulator
// pre-shift); changing one is a breaking change for every consumer.

#include <cstdint>
#include <limits>

#include <gtest/gtest.h>

#include "tap/dsp/sample_traits.h"

namespace {

    using q15 = tap::dsp::sample_traits<std::int16_t>;
    using q31 = tap::dsp::sample_traits<std::int32_t>;
    using f32 = tap::dsp::sample_traits<float>;

    TEST(SampleTraits, CoefficientConversionRoundsAndSaturates) {
        EXPECT_EQ(q15::make_coeff(0.0), 0);
        EXPECT_EQ(q15::make_coeff(1.0), 16384); // Q1.14
        EXPECT_EQ(q15::make_coeff(-1.0), -16384);
        EXPECT_EQ(q15::make_coeff(10.0), 32767); // saturates
        EXPECT_EQ(q15::make_coeff(-10.0), -32768);
        EXPECT_EQ(q31::make_coeff(1.0), 1073741824);   // Q1.30
        EXPECT_EQ(q31::make_coeff(10.0), 2147483647);  // saturates
        EXPECT_FLOAT_EQ(f32::make_coeff(0.25), 0.25f); // unity scale, plain cast
    }

    TEST(SampleTraits, RoundSatRoundsHalfAwayFromZero) {
        using tap::dsp::detail::round_sat;
        EXPECT_EQ(round_sat<std::int16_t>(0.5), 1);
        EXPECT_EQ(round_sat<std::int16_t>(-0.5), -1);
        EXPECT_EQ(round_sat<std::int16_t>(0.49), 0);
        EXPECT_EQ(round_sat<std::int16_t>(-0.49), 0);
        EXPECT_EQ(round_sat<std::int16_t>(1e9), std::numeric_limits<std::int16_t>::max());
        EXPECT_EQ(round_sat<std::int16_t>(-1e9), std::numeric_limits<std::int16_t>::min());
    }

    TEST(SampleTraits, FinalizeSaturates) {
        // Far beyond full scale in the accumulator domain must clamp, not wrap.
        EXPECT_EQ(q15::finalize(std::int64_t{1} << 40), 32767);
        EXPECT_EQ(q15::finalize(-(std::int64_t{1} << 40)), -32768);
        EXPECT_EQ(q31::finalize(std::int64_t{1} << 60), 2147483647);
        EXPECT_EQ(q31::finalize(-(std::int64_t{1} << 60)), -2147483648LL);
    }

    TEST(SampleTraits, FinalizeRoundsHalfUp) {
        // Q29 -> Q15: the rounding constant is 2^13; exactly half rounds up,
        // one below half rounds down, and the same holds on the negative side
        // (round-half-up, not half-away: -0.5 LSB lands on 0).
        EXPECT_EQ(q15::finalize((std::int64_t{1} << 13)), 1);
        EXPECT_EQ(q15::finalize((std::int64_t{1} << 13) - 1), 0);
        EXPECT_EQ(q15::finalize(-(std::int64_t{1} << 13)), 0);
        EXPECT_EQ(q15::finalize(-(std::int64_t{1} << 13) - 1), -1);
        // Q45 -> Q31 uses the identical constant.
        EXPECT_EQ(q31::finalize((std::int64_t{1} << 13)), 1);
        EXPECT_EQ(q31::finalize((std::int64_t{1} << 13) - 1), 0);
    }

    TEST(SampleTraits, Q15MacIsExactInInt64) {
        // Q0.15 x Q1.14 products are exact in int32 and summed without
        // intermediate rounding: full-scale sample times unity coefficient.
        const auto acc = q15::mac(0, std::int16_t{32767}, q15::make_coeff(1.0));
        EXPECT_EQ(acc, std::int64_t{32767} * 16384);
        // finalize returns the sample: one rounding, at the end.
        EXPECT_EQ(q15::finalize(acc), 32767);
    }

    TEST(SampleTraits, Q31MacPreShiftsSixteenBits) {
        // Q0.31 x Q1.30 = Q61 would overflow int64 after ~48 accumulations, so
        // each product drops 16 bits before the add (Q45). Contract: exactly
        // 16, no more (precision), no fewer (headroom).
        const auto acc = q31::mac(0, std::int32_t{1} << 30, std::int32_t{1} << 30);
        EXPECT_EQ(acc, std::int64_t{1} << 44); // (2^60) >> 16
        // Headroom check: ~80 full-scale accumulations stay inside int64.
        std::int64_t big = 0;
        for (int i = 0; i < 80; ++i) {
            big = q31::mac(big, std::numeric_limits<std::int32_t>::max(), std::numeric_limits<std::int32_t>::max());
        }
        EXPECT_GT(big, 0); // no wrap
        EXPECT_EQ(q31::finalize(big), std::numeric_limits<std::int32_t>::max());
    }

    TEST(SampleTraits, FloatMacAccumulatesInDouble) {
        static_assert(std::is_same_v<f32::accum, double>);
        // A double accumulator must hold a contribution a float one would drop:
        // 1.0 + 2^-30 survives in double, vanishes in float.
        const double acc = f32::mac(f32::mac(0.0, 1.0f, 1.0f), 0x1p-15f, 0x1p-15f);
        EXPECT_GT(acc, 1.0);
        EXPECT_DOUBLE_EQ(acc, 1.0 + 0x1p-30);
    }

    TEST(SampleTraits, SilenceIsZero) {
        EXPECT_EQ(q15::silence(), 0);
        EXPECT_EQ(q31::silence(), 0);
        EXPECT_EQ(f32::silence(), 0.0f);
    }

    // The core concept is the kernels' requirement set; all three shipping
    // formats satisfy it (also statically asserted in the header).
    static_assert(tap::dsp::sample_type<float>);
    static_assert(tap::dsp::sample_type<std::int16_t>);
    static_assert(tap::dsp::sample_type<std::int32_t>);
    static_assert(!tap::dsp::sample_type<double>); // no specialization on purpose:
    // double is the golden-model/accumulator domain, not an I/O sample format.

} // namespace
