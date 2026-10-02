// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Contract battery for the sample-format traits, ported and extended from
// SampleRateTap's test_fixed_point.cpp. Every number here is a documented
// contract point (Q formats, rounding mode, saturation, accumulator
// pre-shift); changing one is a breaking change for every consumer.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <type_traits>
#include <utility>

#include <gtest/gtest.h>

#include "tap/dsp/sample_traits.h"

namespace {

    using q15 = tap::dsp::sample_traits<std::int16_t>;
    using q31 = tap::dsp::sample_traits<std::int32_t>;
    using f32 = tap::dsp::sample_traits<float>;
    using f64 = tap::dsp::sample_traits<double>;

    TEST(SampleTraits, CoefficientConversionRoundsAndSaturates) {
        EXPECT_EQ(q15::make_coeff(0.0), 0);
        EXPECT_EQ(q15::make_coeff(1.0), 16384); // Q1.14
        EXPECT_EQ(q15::make_coeff(-1.0), -16384);
        EXPECT_EQ(q15::make_coeff(10.0), 32767); // saturates
        EXPECT_EQ(q15::make_coeff(-10.0), -32768);
        EXPECT_EQ(q31::make_coeff(1.0), 1073741824);   // Q1.30
        EXPECT_EQ(q31::make_coeff(10.0), 2147483647);  // saturates
        EXPECT_FLOAT_EQ(f32::make_coeff(0.25), 0.25f); // unity scale, plain cast
        EXPECT_EQ(f64::make_coeff(0.1), 0.1);          // identity
    }

    TEST(SampleTraits, CoefficientScaleIsTheFormatsUnity) {
        // k_coeff_scale is derived from the named fraction-bit constants and
        // must be the number make_coeff maps 1.0 to; quantize.h multiplies by
        // one and compares against the other.
        static_assert(q15::k_coeff_frac_bits == 14);
        static_assert(q31::k_coeff_frac_bits == 30);
        EXPECT_EQ(q15::k_coeff_scale, 16384.0);
        EXPECT_EQ(q31::k_coeff_scale, 1073741824.0);
        EXPECT_EQ(static_cast<double>(q15::make_coeff(1.0)), q15::k_coeff_scale);
        EXPECT_EQ(static_cast<double>(q31::make_coeff(1.0)), q31::k_coeff_scale);
        EXPECT_EQ(f32::k_coeff_scale, 1.0);
        EXPECT_EQ(f64::k_coeff_scale, 1.0);
    }

    TEST(SampleTraits, QLadderClosesFromTheNamedConstants) {
        // The shifts are derived, not written: accumulator format = sample +
        // coefficient fraction bits - pre-shift; finalize = accumulator ->
        // sample. Both profiles land on the documented 14-bit rounding.
        static_assert(q15::k_sample_frac_bits == 15 && q15::k_accum_pre_shift == 0);
        static_assert(q15::k_accum_frac_bits == 29 && q15::k_finalize_shift == 14);
        static_assert(q31::k_sample_frac_bits == 31 && q31::k_accum_pre_shift == 16);
        static_assert(q31::k_accum_frac_bits == 45 && q31::k_finalize_shift == 14);
        static_assert(q15::k_is_fixed_point && q31::k_is_fixed_point);
        static_assert(!f32::k_is_fixed_point && !f64::k_is_fixed_point);
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
        // Q45 -> Q31 uses the identical constant, on both sides of zero.
        EXPECT_EQ(q31::finalize((std::int64_t{1} << 13)), 1);
        EXPECT_EQ(q31::finalize((std::int64_t{1} << 13) - 1), 0);
        EXPECT_EQ(q31::finalize(-(std::int64_t{1} << 13)), 0);
        EXPECT_EQ(q31::finalize(-(std::int64_t{1} << 13) - 1), -1);
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

    TEST(SampleTraits, DoubleIsTheIdentityDatapath) {
        static_assert(std::is_same_v<f64::coeff, double>);
        static_assert(std::is_same_v<f64::accum, double>);
        // Nothing rounds, nothing saturates: the golden model passes values
        // through untouched, so a traits-based primitive's double profile is
        // the same arithmetic as its hand-written double sibling.
        EXPECT_EQ(f64::make_coeff(1e300), 1e300);
        EXPECT_EQ(f64::finalize(1e300), 1e300);
        EXPECT_EQ(f64::mac(1.0, 0x1p-40, 1.0), 1.0 + 0x1p-40);
        EXPECT_EQ(f64::silence(), 0.0);
    }

    TEST(SampleTraits, SilenceIsZero) {
        EXPECT_EQ(q15::silence(), 0);
        EXPECT_EQ(q31::silence(), 0);
        EXPECT_EQ(f32::silence(), 0.0f);
        EXPECT_EQ(f64::silence(), 0.0);
    }

    // The core concept is the kernels' requirement set; all four shipping
    // formats satisfy it (also statically asserted in the header). double is
    // a sample format because it is the golden model of every primitive
    // (CLAUDE.md): a traits-based primitive instantiates its reference
    // profile through the same substrate its embedded profiles use, and the
    // cross-precision tests measure Q15/Q31/float against it.
    static_assert(tap::dsp::sample_type<double>);
    static_assert(tap::dsp::sample_type<float>);
    static_assert(tap::dsp::sample_type<std::int16_t>);
    static_assert(tap::dsp::sample_type<std::int32_t>);
    static_assert(!tap::dsp::sample_type<std::int8_t>); // no specialization: not a format
    static_assert(!tap::dsp::sample_type<long double>);

    // Every operation is constexpr: the concept's additive-identity
    // requirement is checked at compile time in the header, and here for the
    // value the kernels start from.
    static_assert(q15::mac(q15::accum{}, std::int16_t{-32768}, std::int16_t{-16384}) == std::int64_t{32768} * 16384);
    static_assert(f64::mac(f64::accum{}, 3.0, -2.0) == -6.0);

    // ------------------------------------------------------------------
    // finalize_divided<S, D>: finalize() of acc / D under one rounding.

    /// The exact round-half-up quotient floor((acc + D 2^(f-1)) / (D 2^f)).
    std::int64_t exact_quotient(std::int64_t acc, std::int64_t divisor, int f) {
        const std::int64_t den = divisor << f;
        const std::int64_t num = acc + den / 2;
        std::int64_t       q   = num / den;
        if (num % den != 0 && num < 0) {
            --q;
        }
        return q;
    }

    template <std::uint32_t D>
    void expect_q15_multiples_exact() {
        for (std::int32_t x = -32768; x <= 32767; ++x) {
            const std::int64_t acc = static_cast<std::int64_t>(D) * (std::int64_t{1} << q15::k_finalize_shift) * x;
            ASSERT_EQ((tap::dsp::finalize_divided<std::int16_t, D>(acc)), x) << "D=" << D;
        }
    }

    template <std::uint32_t D>
    void expect_q31_multiples_exact() {
        for (std::int64_t x : {std::int64_t{-2147483648LL}, std::int64_t{-1}, std::int64_t{0}, std::int64_t{1},
                               std::int64_t{123456789}, std::int64_t{2147483647}}) {
            const std::int64_t acc = static_cast<std::int64_t>(D) * (std::int64_t{1} << q31::k_finalize_shift) * x;
            EXPECT_EQ((tap::dsp::finalize_divided<std::int32_t, D>(acc)), x) << "D=" << D;
        }
    }

    // Every exact multiple of D 2^f — the accumulator of DC through rows
    // summing to D x unity — comes out exact: all 65536 Q15 values for each
    // divisor of the family's vocabulary and two odd ones, a sample of Q31.
    TEST(SampleTraits, FinalizeDividedIsExactAtEveryMultipleOfTheDivisor) {
        expect_q15_multiples_exact<1>();
        expect_q15_multiples_exact<2>();
        expect_q15_multiples_exact<3>();
        expect_q15_multiples_exact<4>();
        expect_q15_multiples_exact<5>();
        expect_q15_multiples_exact<6>();
        expect_q15_multiples_exact<7>();
        expect_q15_multiples_exact<8>();
        expect_q31_multiples_exact<1>();
        expect_q31_multiples_exact<3>();
        expect_q31_multiples_exact<6>();
        expect_q31_multiples_exact<8>();
    }

    // D = 1 is finalize(); a power of two is finalize's round-half-up with the
    // shift widened (exact), Q31's other divisors the exact quotient, and
    // Q15's multiply-back the exact quotient except within D 2^-11 of a
    // rounding boundary, where it may land one step away.
    template <std::uint32_t D>
    void expect_q15_quotient(std::uint64_t seed) {
        constexpr int f = q15::k_finalize_shift;
        for (int i = 0; i < 20000; ++i) {
            seed                     = seed * 6364136223846793005ULL + 1442695040888963407ULL;
            const auto         acc   = static_cast<std::int64_t>(seed >> 30) - (std::int64_t{1} << 33);
            const std::int64_t exact = exact_quotient(acc, D, f);
            const auto         got   = static_cast<std::int64_t>(tap::dsp::finalize_divided<std::int16_t, D>(acc));
            const std::int64_t clamp = std::min<std::int64_t>(32767, std::max<std::int64_t>(-32768, exact));
            if constexpr ((D & (D - 1)) == 0) {
                ASSERT_EQ(got, clamp) << "D=" << D << " acc=" << acc;
            }
            else {
                // Distance of the true quotient from its rounding boundary.
                const double quotient = static_cast<double>(acc) / (static_cast<double>(D) * (1 << f));
                const double frac     = quotient - std::floor(quotient);
                if (std::abs(frac - 0.5) > static_cast<double>(D) / 2048.0) {
                    ASSERT_EQ(got, clamp) << "D=" << D << " acc=" << acc;
                }
                else {
                    ASSERT_LE(std::abs(got - clamp), 1) << "D=" << D << " acc=" << acc;
                }
            }
        }
    }

    TEST(SampleTraits, FinalizeDividedRoundsTheQuotientHalfUp) {
        for (std::int64_t acc : {std::int64_t{0}, std::int64_t{8191}, std::int64_t{8192}, std::int64_t{-8192},
                                 std::int64_t{-8193}, std::int64_t{1} << 40, -(std::int64_t{1} << 40)}) {
            EXPECT_EQ((tap::dsp::finalize_divided<std::int16_t, 1>(acc)), q15::finalize(acc)) << acc;
            EXPECT_EQ((tap::dsp::finalize_divided<std::int32_t, 1>(acc)), q31::finalize(acc)) << acc;
            EXPECT_EQ((tap::dsp::finalize_divided<double, 1>(static_cast<double>(acc))),
                      f64::finalize(static_cast<double>(acc)));
        }
        expect_q15_quotient<2>(1);
        expect_q15_quotient<3>(2);
        expect_q15_quotient<6>(3);
        expect_q15_quotient<8>(4);
        std::uint64_t seed = 5;
        for (int i = 0; i < 20000; ++i) {
            seed           = seed * 6364136223846793005ULL + 1442695040888963407ULL;
            const auto acc = static_cast<std::int64_t>(seed >> 14) - (std::int64_t{1} << 49);
            for (const auto& [got, d] :
                 {std::pair{static_cast<std::int64_t>(tap::dsp::finalize_divided<std::int32_t, 3>(acc)), 3},
                  std::pair{static_cast<std::int64_t>(tap::dsp::finalize_divided<std::int32_t, 6>(acc)), 6},
                  std::pair{static_cast<std::int64_t>(tap::dsp::finalize_divided<std::int32_t, 8>(acc)), 8}}) {
                const std::int64_t exact = exact_quotient(acc, d, q31::k_finalize_shift);
                ASSERT_EQ(got, std::min<std::int64_t>(2147483647, std::max<std::int64_t>(-2147483648LL, exact)))
                    << "D=" << d << " acc=" << acc;
            }
        }
        // Floating: finalize of the quotient in the accumulator domain.
        EXPECT_EQ((tap::dsp::finalize_divided<float, 3>(1.0)), static_cast<float>(1.0 / 3.0));
        EXPECT_EQ((tap::dsp::finalize_divided<double, 6>(-3.0)), -0.5);
    }

    TEST(SampleTraits, FinalizeDividedSaturates) {
        const std::int64_t big = std::int64_t{1} << 36;
        EXPECT_EQ((tap::dsp::finalize_divided<std::int16_t, 3>(big)), 32767);
        EXPECT_EQ((tap::dsp::finalize_divided<std::int16_t, 3>(-big)), -32768);
        EXPECT_EQ((tap::dsp::finalize_divided<std::int16_t, 8>(big)), 32767);
        EXPECT_EQ((tap::dsp::finalize_divided<std::int32_t, 6>(std::int64_t{1} << 62)), 2147483647);
        EXPECT_EQ((tap::dsp::finalize_divided<std::int32_t, 6>(-(std::int64_t{1} << 62))), -2147483648LL);
    }

} // namespace
