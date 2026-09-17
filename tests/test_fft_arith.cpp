// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Contract battery for the butterfly arithmetic trait (fft/fft_arith.h),
// the specification the fixed-point FFT kernel is written against. Every
// number here is a documented contract point: the Q1.30 twiddle format, the
// single rounding point of the multiply, saturation at the int32 rails, the
// Q15 profile's two guard bits, and the block-headroom definition.

#include <array>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/fft/fft_arith.h"

namespace {

    using q31 = tap::dsp::fft_arith<std::int32_t>;
    using q15 = tap::dsp::fft_arith<std::int16_t>;
    using f32 = tap::dsp::fft_arith<float>;
    using f64 = tap::dsp::fft_arith<double>;

    constexpr std::int32_t k_i32_max = std::numeric_limits<std::int32_t>::max();
    constexpr std::int32_t k_i32_min = std::numeric_limits<std::int32_t>::min();
    constexpr std::int16_t k_i16_max = std::numeric_limits<std::int16_t>::max();
    constexpr std::int16_t k_i16_min = std::numeric_limits<std::int16_t>::min();

    TEST(FftArith, TwiddleUnityIsRepresentable) {
        // Q1.30 for both fixed profiles, one type, 1.0 = 2^30 exactly.
        static_assert(std::is_same_v<q31::coeff, std::int32_t>);
        static_assert(std::is_same_v<q15::coeff, q31::coeff>);
        static_assert(q31::k_coeff_frac_bits == 30 && q15::k_coeff_frac_bits == 30);
        EXPECT_EQ(q31::make_coeff(1.0), 1073741824);
        EXPECT_EQ(q15::make_coeff(1.0), 1073741824);
        EXPECT_EQ(q31::make_coeff(-1.0), -1073741824);
        EXPECT_EQ(q31::make_coeff(0.5), 536870912);
        EXPECT_EQ(q31::k_coeff_one, 1073741824);
        // Rounded once, half away from zero, saturating (round_sat).
        EXPECT_EQ(q31::make_coeff(0x1p-31), 1); // 0.5 LSB rounds up
        EXPECT_EQ(q31::make_coeff(-0x1p-31), -1);
        EXPECT_EQ(q31::make_coeff(3.0), k_i32_max);
        // A unity twiddle is the identity on any sample.
        EXPECT_EQ(q31::mul_coeff(123456789, q31::k_coeff_one), 123456789);
        EXPECT_EQ(q31::mul_coeff(k_i32_max, q31::k_coeff_one), k_i32_max);
        EXPECT_EQ(q31::mul_coeff(k_i32_min, q31::k_coeff_one), k_i32_min);
    }

    TEST(FftArith, MulCoeffIsExactOnRepresentableProducts) {
        // x * w in Qm.n x Q1.30: when the exact product has no bits below the
        // shift, the result is the exact product (hand-computed).
        EXPECT_EQ(q31::mul_coeff(1 << 20, 1 << 29), 1 << 19); // 2^20 * 0.5
        EXPECT_EQ(q31::mul_coeff(-(1 << 20), 1 << 29), -(1 << 19));
        EXPECT_EQ(q31::mul_coeff(1 << 20, -(1 << 29)), -(1 << 19));
        EXPECT_EQ(q31::mul_coeff(3 << 10, 1 << 28), 3 << 8); // 3072 * 0.25 = 768
        EXPECT_EQ(q31::mul_coeff(1000, 0), 0);
        EXPECT_EQ(q31::mul_coeff(0, 1 << 30), 0);
        EXPECT_EQ(q31::mul_coeff(1 << 30, 1 << 30), 1 << 30); // 0.5 * 1.0 in Q0.31
        EXPECT_EQ(q31::mul_coeff(1 << 30, 3 << 28), 3 << 28); // 0.5 * 0.75 = 0.375
        // The float and double profiles are the plain product.
        EXPECT_EQ(f32::mul_coeff(0.5f, 0.75f), 0.375f);
        EXPECT_EQ(f64::mul_coeff(0.5, 0.75), 0.375);
    }

    TEST(FftArith, MulCoeffRoundsHalfUp) {
        // 1 * 0.5 in Q1.30 = 2^29 in the product: exactly half an LSB, which
        // rounds up to 1 for the positive product and to 0 for the negative
        // product (half-up, not half-away).
        EXPECT_EQ(q31::mul_coeff(1, 1 << 29), 1);
        EXPECT_EQ(q31::mul_coeff(-1, 1 << 29), 0);
        EXPECT_EQ(q31::mul_coeff(1, -(1 << 29)), 0);
        // Just below half rounds down; just above rounds up, on both sides.
        EXPECT_EQ(q31::mul_coeff(1, (1 << 29) - 1), 0);
        EXPECT_EQ(q31::mul_coeff(1, (1 << 29) + 1), 1);
        EXPECT_EQ(q31::mul_coeff(-1, (1 << 29) - 1), 0);
        EXPECT_EQ(q31::mul_coeff(-1, (1 << 29) + 1), -1);
        // 3 * 0.5 = 1.5 -> 2; -3 * 0.5 = -1.5 -> -1.
        EXPECT_EQ(q31::mul_coeff(3, 1 << 29), 2);
        EXPECT_EQ(q31::mul_coeff(-3, 1 << 29), -1);
        // One rounding, not two: 2^30 - 1 (0.999999999) times 2^30 - 1.
        // Exact product = (2^30 - 1)^2 = 2^60 - 2^31 + 1; >> 30 with half-up:
        // floor((2^60 - 2^31 + 1 + 2^29) / 2^30) = 2^30 - 2 + floor((1 + 2^29) / 2^30) = 2^30 - 2.
        EXPECT_EQ(q31::mul_coeff((1 << 30) - 1, (1 << 30) - 1), (1 << 30) - 2);
    }

    TEST(FftArith, MulCoeffSaturatesAtTheRail) {
        // With |w| <= 1.0 the only saturating input is INT32_MIN * -1.0, whose
        // exact product is +2^31.
        EXPECT_EQ(q31::mul_coeff(k_i32_min, -q31::k_coeff_one), k_i32_max);
        EXPECT_EQ(q31::mul_coeff(k_i32_min, -q31::k_coeff_one + 1), k_i32_max - 1); // exact 2^31 - 2 + 2^-30: fits
        EXPECT_EQ(q31::mul_coeff(k_i32_max, -q31::k_coeff_one), k_i32_min + 1);
        EXPECT_EQ(q31::mul_coeff(k_i32_min, q31::k_coeff_one), k_i32_min);
        // Out-of-contract twiddles (|w| > 1.0) saturate rather than wrap.
        EXPECT_EQ(q31::mul_coeff(k_i32_max, k_i32_max), k_i32_max);
        EXPECT_EQ(q31::mul_coeff(k_i32_min, k_i32_max), k_i32_min);
        EXPECT_EQ(q31::mul_coeff(k_i32_min, k_i32_min), k_i32_max);
    }

    TEST(FftArith, AddSubSaturateAtTheRails) {
        EXPECT_EQ(q31::add(1, 2), 3);
        EXPECT_EQ(q31::sub(1, 2), -1);
        EXPECT_EQ(q31::add(k_i32_max, 1), k_i32_max);
        EXPECT_EQ(q31::add(k_i32_max, k_i32_max), k_i32_max);
        EXPECT_EQ(q31::add(k_i32_min, -1), k_i32_min);
        EXPECT_EQ(q31::add(k_i32_min, k_i32_min), k_i32_min);
        EXPECT_EQ(q31::sub(k_i32_min, 1), k_i32_min);
        EXPECT_EQ(q31::sub(k_i32_max, -1), k_i32_max);
        EXPECT_EQ(q31::sub(0, k_i32_min), k_i32_max); // -INT32_MIN does not exist
        EXPECT_EQ(q31::sub(k_i32_min, k_i32_max), k_i32_min);
        // In range, the plain results.
        EXPECT_EQ(q31::add(k_i32_max, k_i32_min), -1);
        EXPECT_EQ(q31::sub(k_i32_min, k_i32_min), 0);
        EXPECT_EQ(f64::add(0.5, 0.25), 0.75);
        EXPECT_EQ(f32::sub(0.5f, 0.25f), 0.25f);
    }

    TEST(FftArith, ShrRoundRoundsHalfUp) {
        EXPECT_EQ(q31::shr_round(8, 2), 2);
        EXPECT_EQ(q31::shr_round(2, 2), 1);   // 0.5 -> 1
        EXPECT_EQ(q31::shr_round(-2, 2), 0);  // -0.5 -> 0 (half-up)
        EXPECT_EQ(q31::shr_round(1, 2), 0);   // 0.25 -> 0
        EXPECT_EQ(q31::shr_round(-1, 2), 0);  // -0.25 -> 0
        EXPECT_EQ(q31::shr_round(3, 2), 1);   // 0.75 -> 1
        EXPECT_EQ(q31::shr_round(-3, 2), -1); // -0.75 -> -1
        EXPECT_EQ(q31::shr_round(-6, 2), -1); // -1.5 -> -1
        EXPECT_EQ(q31::shr_round(6, 2), 2);   // 1.5 -> 2
        EXPECT_EQ(q31::shr_round(7, 0), 7);   // identity
        // The rounding add cannot overflow: full scale by one bit.
        EXPECT_EQ(q31::shr_round(k_i32_max, 1), 1 << 30); // (2^31 - 1 + 1) / 2
        EXPECT_EQ(q31::shr_round(k_i32_min, 1), -(1 << 30));
        EXPECT_EQ(q31::shr_round(k_i32_min, 31), -1);
        EXPECT_EQ(q31::shr_round(k_i32_max, 31), 1);      // 0.99999 -> 1
        EXPECT_EQ(q31::shr_round(k_i32_min + 1, 31), -1); // -0.99999 -> -1
        // Floating: exact power-of-two scale.
        EXPECT_EQ(f64::shr_round(3.0, 2), 0.75);
        EXPECT_EQ(f32::shr_round(-3.0f, 1), -1.5f);
        EXPECT_EQ(f64::shr_round(1.0, 0), 1.0);
    }

    TEST(FftArith, HeadroomBitsOnZeroOneAndFullScale) {
        const std::array<std::int32_t, 4> zeros{0, 0, 0, 0};
        EXPECT_EQ(q31::headroom_bits(zeros.data(), zeros.size()), 31);
        EXPECT_EQ(q31::headroom_bits(zeros.data(), 0), 31); // empty block
        const std::array<std::int32_t, 1> one{1};
        EXPECT_EQ(q31::headroom_bits(one.data(), one.size()), 30);
        const std::array<std::int32_t, 1> minus_one{-1};
        EXPECT_EQ(q31::headroom_bits(minus_one.data(), 1), 31); // -1 << 31 still fits
        const std::array<std::int32_t, 1> minus_two{-2};
        EXPECT_EQ(q31::headroom_bits(minus_two.data(), 1), 30);
        const std::array<std::int32_t, 3> full{0, k_i32_max, 0};
        EXPECT_EQ(q31::headroom_bits(full.data(), full.size()), 0);
        const std::array<std::int32_t, 3> full_neg{0, k_i32_min, 0};
        EXPECT_EQ(q31::headroom_bits(full_neg.data(), full_neg.size()), 0);
        const std::array<std::int32_t, 2> half{1 << 30, -(1 << 30)};
        EXPECT_EQ(q31::headroom_bits(half.data(), half.size()), 0); // 2^30 << 1 = 2^31 does not fit
        const std::array<std::int32_t, 2> quarter{(1 << 29) + 5, -(1 << 29)};
        EXPECT_EQ(q31::headroom_bits(quarter.data(), quarter.size()), 1);
        // The definition, checked on the whole bit ladder: x = 1 << (30 - s)
        // has headroom s, and -x has s + 1 (two's complement holds one more
        // negative power of two: -2^k << (31 - k) is INT32_MIN and fits).
        for (int s = 0; s <= 30; ++s) {
            const std::int32_t x = std::int32_t{1} << (30 - s);
            EXPECT_EQ(q31::headroom_bits(&x, 1), s) << "s=" << s;
            const std::int32_t nx = -x;
            EXPECT_EQ(q31::headroom_bits(&nx, 1), s + 1) << "s=" << s;
        }
        // Floating blocks report no headroom: no block scaling.
        const std::array<double, 2> fd{0.5, 1e-9};
        EXPECT_EQ(f64::headroom_bits(fd.data(), fd.size()), 0);
    }

    TEST(FftArith, WidenPlacesTwoGuardBits) {
        static_assert(q15::k_guard_bits == 2 && q15::k_widen_shift == 14);
        static_assert(std::is_same_v<q15::wide, std::int32_t>);
        static_assert(std::is_same_v<q15::work, q31>);
        EXPECT_EQ(q15::widen(1), 1 << 14);
        EXPECT_EQ(q15::widen(k_i16_max), 32767 << 14);
        EXPECT_EQ(q15::widen(k_i16_min), -(1 << 29));
        // Full scale sits two bits below the int32 sign: headroom exactly 2.
        const std::int32_t fs = q15::widen(k_i16_max);
        EXPECT_EQ(q31::headroom_bits(&fs, 1), 2);
        const std::int32_t nfs = q15::widen(k_i16_min);
        EXPECT_EQ(q31::headroom_bits(&nfs, 1), 2);
        // 4x growth from full scale fits (the fixed-scaling argument).
        EXPECT_EQ(q31::add(q31::add(nfs, nfs), q31::add(nfs, nfs)), k_i32_min);
        EXPECT_EQ(q31::add(q31::add(fs, fs), q31::add(fs, fs)), 4 * (32767 << 14));
    }

    TEST(FftArith, NarrowRoundsHalfUpAndSaturates) {
        EXPECT_EQ(q15::narrow(1 << 14), 1);
        EXPECT_EQ(q15::narrow(1 << 13), 1); // exactly half rounds up
        EXPECT_EQ(q15::narrow((1 << 13) - 1), 0);
        EXPECT_EQ(q15::narrow(-(1 << 13)), 0); // -0.5 LSB -> 0 (half-up)
        EXPECT_EQ(q15::narrow(-(1 << 13) - 1), -1);
        EXPECT_EQ(q15::narrow(k_i32_max), k_i16_max);
        EXPECT_EQ(q15::narrow(k_i32_min), k_i16_min);
        EXPECT_EQ(q15::narrow(32768 << 14), k_i16_max); // one LSB over full scale saturates
        EXPECT_EQ(q15::narrow((32767 << 14) + (1 << 13)), k_i16_max);
        EXPECT_EQ(q15::narrow(-(32769 << 14)), k_i16_min);
        EXPECT_EQ(q15::narrow(-(1 << 29)), k_i16_min);
    }

    TEST(FftArith, WidenNarrowRoundTripsEveryInt16) {
        for (int v = k_i16_min; v <= k_i16_max; ++v) {
            const auto x = static_cast<std::int16_t>(v);
            ASSERT_EQ(q15::narrow(q15::widen(x)), x) << v;
        }
        // And the identities of the other profiles.
        EXPECT_EQ(q31::narrow(q31::widen(k_i32_min)), k_i32_min);
        EXPECT_EQ(f32::narrow(f32::widen(0.5f)), 0.5f);
        EXPECT_EQ(f64::narrow(f64::widen(-0.25)), -0.25);
    }

    TEST(FftArith, EveryOperationIsConstexprAndNoexcept) {
        static_assert(noexcept(q31::mul_coeff(1, 1)) && noexcept(q31::add(1, 1)) && noexcept(q31::sub(1, 1)));
        static_assert(noexcept(q31::shr_round(1, 1)) && noexcept(q31::headroom_bits(nullptr, 0)));
        static_assert(noexcept(q15::widen(1)) && noexcept(q15::narrow(1)) && noexcept(q31::make_coeff(1.0)));
        static_assert(q31::mul_coeff(1 << 20, 1 << 29) == 1 << 19);
        static_assert(q31::add(k_i32_max, 1) == k_i32_max);
        static_assert(q31::shr_round(-2, 2) == 0);
        static_assert(q15::narrow(q15::widen(-12345)) == -12345);
        constexpr std::array<std::int32_t, 2> block{1 << 20, -(1 << 24)};
        static_assert(q31::headroom_bits(block.data(), block.size()) == 7); // -2^24 << 7 = INT32_MIN
        static_assert(f64::shr_round(1.0, 3) == 0.125);
    }

    // A butterfly on the widened Q15 profile, hand-computed end to end: the
    // arithmetic the kernel will compose, checked once here so the pieces
    // are known to fit together (formats, rounding point, narrowing).
    TEST(FftArith, WidenedButterflyHandComputed) {
        // x0 = 0.5, x1 = 0.25 in Q0.15; w = cos(60 deg) = 0.5 in Q1.30.
        const std::int32_t a = q15::widen(16384);
        const std::int32_t b = q15::widen(8192);
        const auto         w = q15::make_coeff(0.5);
        const std::int32_t t = q31::mul_coeff(b, w); // 0.125 in Q2.29
        EXPECT_EQ(t, 4096 << 14);
        const std::int32_t y0 = q31::add(a, t); // 0.625
        const std::int32_t y1 = q31::sub(a, t); // 0.375
        EXPECT_EQ(q15::narrow(y0), 20480);
        EXPECT_EQ(q15::narrow(y1), 12288);
        // Scaled by one bit before narrowing (fixed scaling): half those.
        EXPECT_EQ(q15::narrow(q31::shr_round(y0, 1)), 10240);
        EXPECT_EQ(q15::narrow(q31::shr_round(y1, 1)), 6144);
    }

} // namespace
