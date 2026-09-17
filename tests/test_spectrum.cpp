// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Locks down the tap::dsp::packed_spectrum view. Against real transforms
// (float/double): the DC and Nyquist slots, bin k at a[2k] / a[2k+1], the sign
// of im(k) for an on-bin sine under the W = exp(+2*pi*i/N) convention, power()
// at the edges and in the interior, Parseval over the packing, the
// engineering-convention accessor's conjugation, and that a value written
// through the view inverts to the documented tone. Over hand-filled buffers
// (float/double/int16/int32): the slots, write-through, deduction, and
// power()'s promotion before the multiply.

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <random>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "tap/dsp/fft.h"
#include "tap/dsp/fft/spectrum.h"

namespace {

    template <typename Sample>
    constexpr double k_tolerance = 0.0;
    template <>
    constexpr double k_tolerance<double> = 1e-12;
    template <>
    constexpr double k_tolerance<float> = 2e-5;

    template <typename Sample>
    std::vector<Sample> forward(std::vector<Sample> x) {
        tap::dsp::basic_real_fft<Sample> fft(x.size());
        fft.forward_inplace(x.data());
        return x;
    }

    template <typename Sample>
    std::vector<Sample> tone(std::size_t n, std::size_t bin, bool sine) {
        std::vector<Sample> x(n);
        const double        w = 2.0 * std::numbers::pi * static_cast<double>(bin) / static_cast<double>(n);
        for (std::size_t j = 0; j < n; ++j) {
            const double arg = w * static_cast<double>(j);
            x[j]             = static_cast<Sample>(sine ? std::sin(arg) : std::cos(arg));
        }
        return x;
    }

    template <typename Sample>
    std::vector<Sample> random_signal(std::size_t n, unsigned seed) {
        std::mt19937                           gen(seed);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        std::vector<Sample>                    x(n);
        for (auto& v : x) {
            v = static_cast<Sample>(dist(gen));
        }
        return x;
    }

    // Transform-backed battery: the two profiles basic_real_fft ships today.
    template <typename Sample>
    class packed_spectrum_test : public ::testing::Test {};

    using sample_types = ::testing::Types<float, double>;
    TYPED_TEST_SUITE(packed_spectrum_test, sample_types);

    // Hand-filled battery: every value type the view admits, the fixed-point
    // profiles included, so Stage 3 does not have to reopen this header.
    template <typename Sample>
    class packed_spectrum_slot_test : public ::testing::Test {};

    using slot_types = ::testing::Types<float, double, std::int16_t, std::int32_t>;
    TYPED_TEST_SUITE(packed_spectrum_slot_test, slot_types);

    TYPED_TEST(packed_spectrum_slot_test, SizeAndBinCount) {
        std::vector<TypeParam>                           a(1024, TypeParam(0));
        const tap::dsp::packed_spectrum<const TypeParam> s(a.data(), a.size());
        EXPECT_EQ(s.size(), 1024u);
        EXPECT_EQ(s.num_bins(), 513u);
        EXPECT_EQ(s.data(), a.data());
    }

    TYPED_TEST(packed_spectrum_test, DcConstantLandsInDc) {
        constexpr std::size_t n   = 64;
        const auto            a   = forward(std::vector<TypeParam>(n, TypeParam(1)));
        const double          tol = k_tolerance<TypeParam> * n;

        const tap::dsp::packed_spectrum<const TypeParam> s(a.data(), n);
        EXPECT_NEAR(s.dc(), static_cast<double>(n), tol);
        EXPECT_NEAR(s.nyquist(), 0.0, tol);
        for (std::size_t k = 1; k < n / 2; ++k) {
            EXPECT_NEAR(s.re(k), 0.0, tol) << "bin " << k;
            EXPECT_NEAR(s.im(k), 0.0, tol) << "bin " << k;
        }
        EXPECT_NEAR(s.power(0), static_cast<double>(n) * n, tol * n);
    }

    TYPED_TEST(packed_spectrum_test, AlternatingSignLandsInNyquist) {
        constexpr std::size_t  n = 64;
        std::vector<TypeParam> x(n);
        for (std::size_t j = 0; j < n; ++j) {
            x[j] = (j % 2 == 0) ? TypeParam(1) : TypeParam(-1);
        }
        const auto   a   = forward(x);
        const double tol = k_tolerance<TypeParam> * n;

        const tap::dsp::packed_spectrum<const TypeParam> s(a.data(), n);
        EXPECT_NEAR(s.dc(), 0.0, tol);
        EXPECT_NEAR(s.nyquist(), static_cast<double>(n), tol);
        EXPECT_NEAR(s.power(n / 2), static_cast<double>(n) * n, tol * n);
        EXPECT_NEAR(s.power(s.num_bins() - 1), static_cast<double>(n) * n, tol * n);
    }

    // The documented sign convention (W = exp(+2*pi*i/N)): an on-bin cosine
    // is real (+N/2), an on-bin sine is imaginary with a PLUS sign (+N/2) in
    // the native reading, and MINUS in the engineering reading.
    TYPED_TEST(packed_spectrum_test, OnBinCosineAndSinePinTheSignOfIm) {
        constexpr std::size_t n    = 128;
        constexpr std::size_t bin  = 5;
        const auto            c    = forward(tone<TypeParam>(n, bin, false));
        const auto            s    = forward(tone<TypeParam>(n, bin, true));
        const double          half = static_cast<double>(n) / 2.0;
        const double          tol  = k_tolerance<TypeParam> * n;

        const tap::dsp::packed_spectrum<const TypeParam> cosine(c.data(), n);
        const tap::dsp::packed_spectrum<const TypeParam> sine(s.data(), n);
        EXPECT_NEAR(cosine.re(bin), half, tol);
        EXPECT_NEAR(cosine.im(bin), 0.0, tol);
        EXPECT_NEAR(sine.re(bin), 0.0, tol);
        EXPECT_NEAR(sine.im(bin), half, tol); // +N/2, not -N/2

        const std::complex<TypeParam> eng = sine.bin_engineering(bin);
        EXPECT_NEAR(eng.real(), 0.0, tol);
        EXPECT_NEAR(eng.imag(), -half, tol); // conjugated: the engineering DFT's -N/2
        EXPECT_EQ(eng.real(), sine.re(bin));
        EXPECT_EQ(eng.imag(), -sine.im(bin));

        // Both tones carry N^2/4 of power in their bin and none elsewhere.
        EXPECT_NEAR(cosine.power(bin), half * half, tol * n);
        EXPECT_NEAR(sine.power(bin), half * half, tol * n);
        for (std::size_t k = 0; k < cosine.num_bins(); ++k) {
            if (k == bin) {
                continue;
            }
            EXPECT_NEAR(cosine.power(k), 0.0, tol * n) << "cos leak, bin " << k;
            EXPECT_NEAR(sine.power(k), 0.0, tol * n) << "sin leak, bin " << k;
        }
    }

    TYPED_TEST(packed_spectrum_slot_test, AccessorsAreTheDocumentedSlots) {
        // Fill a[i] = i so every slot is distinguishable, and read it back
        // through the view: a[0], a[1], a[2k], a[2k+1].
        constexpr std::size_t  n = 32;
        std::vector<TypeParam> a(n);
        for (std::size_t i = 0; i < n; ++i) {
            a[i] = static_cast<TypeParam>(i);
        }
        const tap::dsp::packed_spectrum<const TypeParam> s(a.data(), n);
        EXPECT_EQ(s.dc(), a[0]);
        EXPECT_EQ(s.nyquist(), a[1]);
        for (std::size_t k = 1; k < n / 2; ++k) {
            EXPECT_EQ(s.re(k), a[2 * k]) << "bin " << k;
            EXPECT_EQ(s.im(k), a[2 * k + 1]) << "bin " << k;
            EXPECT_EQ(s.power(k), a[2 * k] * a[2 * k] + a[2 * k + 1] * a[2 * k + 1]) << "bin " << k;
        }
        EXPECT_EQ(s.power(0), a[0] * a[0]);
        EXPECT_EQ(s.power(n / 2), a[1] * a[1]);
    }

    TYPED_TEST(packed_spectrum_slot_test, MutableViewWritesThrough) {
        constexpr std::size_t  n = 16;
        std::vector<TypeParam> a(n, TypeParam(0));

        const tap::dsp::packed_spectrum<TypeParam> s(a.data(), n);
        s.dc()      = TypeParam(1);
        s.nyquist() = TypeParam(2);
        s.re(3)     = TypeParam(3);
        s.im(3) += TypeParam(4);
        EXPECT_EQ(a[0], TypeParam(1));
        EXPECT_EQ(a[1], TypeParam(2));
        EXPECT_EQ(a[6], TypeParam(3));
        EXPECT_EQ(a[7], TypeParam(4));

        // A mutable view converts to the read-only one over the same buffer.
        const tap::dsp::packed_spectrum<const TypeParam> r = s;
        EXPECT_EQ(r.data(), a.data());
        EXPECT_EQ(r.re(3), TypeParam(3));
        static_assert(std::is_same_v<decltype(r.re(3)), const TypeParam&>);
        static_assert(std::is_same_v<decltype(s.re(3)), TypeParam&>);
    }

    TYPED_TEST(packed_spectrum_slot_test, DeductionFollowsThePointer) {
        std::vector<TypeParam>    a(8, TypeParam(0));
        const auto*               ca = a.data();
        tap::dsp::packed_spectrum m(a.data(), a.size());
        tap::dsp::packed_spectrum c(ca, a.size());
        static_assert(std::is_same_v<decltype(m), tap::dsp::packed_spectrum<TypeParam>>);
        static_assert(std::is_same_v<decltype(c), tap::dsp::packed_spectrum<const TypeParam>>);
        EXPECT_EQ(m.num_bins(), 5u);
        EXPECT_EQ(c.num_bins(), 5u);
    }

    TYPED_TEST(packed_spectrum_slot_test, PowerPromotesBeforeTheMultiply) {
        // power_type is the sample type for the floating profiles (the
        // hand-written consumers' exact expression) and int64 for the
        // fixed-point ones, promoted BEFORE the multiply: an int16 pair of
        // 30000 squares to 1.8e9, which int16 arithmetic narrowed to -11776,
        // and an int32 pair of 2^30 squares to 2^61, past int32 entirely.
        using view = tap::dsp::packed_spectrum<const TypeParam>;
        if constexpr (std::is_integral_v<TypeParam>) {
            static_assert(std::is_same_v<typename view::power_type, std::int64_t>);
        }
        else {
            static_assert(std::is_same_v<typename view::power_type, TypeParam>);
        }

        constexpr std::size_t  n = 8;
        std::vector<TypeParam> a(n, TypeParam(0));
        const TypeParam        big = std::is_same_v<TypeParam, std::int16_t>   ? TypeParam(30000)
                                     : std::is_same_v<TypeParam, std::int32_t> ? TypeParam(1 << 30)
                                                                               : TypeParam(30000);
        a[0]                       = big;
        a[1]                       = big;
        a[2]                       = big;
        a[3]                       = big;
        const view s(a.data(), n);

        const std::int64_t b = static_cast<std::int64_t>(big);
        EXPECT_EQ(static_cast<std::int64_t>(s.power(0)), b * b);
        EXPECT_EQ(static_cast<std::int64_t>(s.power(n / 2)), b * b);
        EXPECT_EQ(static_cast<std::int64_t>(s.power(1)), b * b + b * b);
        EXPECT_EQ(s.power(2), typename view::power_type(0));
    }

    TYPED_TEST(packed_spectrum_test, ParsevalHoldsOverThePacking) {
        // The power() docstring's formula, with no hidden factor 2 or 1/N:
        // sum x^2 = (1/N) (power(0) + power(N/2) + 2 sum_{1..N/2-1} power(k)).
        constexpr std::size_t n = 512;
        const auto            x = random_signal<TypeParam>(n, 1234);

        double time_energy = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            time_energy += static_cast<double>(x[i]) * static_cast<double>(x[i]);
        }

        const auto                                       a = forward(x);
        const tap::dsp::packed_spectrum<const TypeParam> s(a.data(), n);
        double freq_energy = static_cast<double>(s.power(0)) + static_cast<double>(s.power(s.num_bins() - 1));
        for (std::size_t k = 1; k < s.num_bins() - 1; ++k) {
            freq_energy += 2.0 * static_cast<double>(s.power(k));
        }
        freq_energy /= static_cast<double>(n);

        const double tol = std::is_same_v<TypeParam, double> ? 1e-9 : 1e-3;
        EXPECT_NEAR(freq_energy, time_energy, tol * time_energy);
    }

    TYPED_TEST(packed_spectrum_test, WritesThroughTheViewInvertToTheDocumentedTones) {
        // The side pvoc writes: set ONLY im(5) = +N/2 through a mutable view
        // and run the real inverse (2/N applied). Under W = exp(+2*pi*i/N)
        // that is a sine, not a negative sine; re(5) = N/2 alone is a cosine.
        constexpr std::size_t n     = 128;
        constexpr std::size_t bin   = 5;
        const double          w     = 2.0 * std::numbers::pi * static_cast<double>(bin) / static_cast<double>(n);
        const double          tol   = k_tolerance<TypeParam> * static_cast<double>(n);
        const TypeParam       scale = TypeParam(2) / static_cast<TypeParam>(n);

        tap::dsp::basic_real_fft<TypeParam> fft(n);

        std::vector<TypeParam> sine(n, TypeParam(0));
        {
            const tap::dsp::packed_spectrum<TypeParam> s(sine.data(), n);
            s.im(bin) = static_cast<TypeParam>(n) / TypeParam(2);
        }
        fft.inverse_inplace(sine.data());

        std::vector<TypeParam> cosine(n, TypeParam(0));
        {
            const tap::dsp::packed_spectrum<TypeParam> c(cosine.data(), n);
            c.re(bin) = static_cast<TypeParam>(n) / TypeParam(2);
        }
        fft.inverse_inplace(cosine.data());

        for (std::size_t j = 0; j < n; ++j) {
            const double arg = w * static_cast<double>(j);
            EXPECT_NEAR(static_cast<double>(sine[j] * scale), std::sin(arg), tol) << "sine, sample " << j;
            EXPECT_NEAR(static_cast<double>(cosine[j] * scale), std::cos(arg), tol) << "cosine, sample " << j;
        }
    }

    TEST(packed_spectrum_contract, AccessorsAreNoexceptAndConstexpr) {
        using view = tap::dsp::packed_spectrum<const float>;
        static_assert(noexcept(std::declval<const view&>().dc()));
        static_assert(noexcept(std::declval<const view&>().nyquist()));
        static_assert(noexcept(std::declval<const view&>().re(1)));
        static_assert(noexcept(std::declval<const view&>().im(1)));
        static_assert(noexcept(std::declval<const view&>().power(1)));
        static_assert(noexcept(std::declval<const view&>().bin_engineering(1)));
        static_assert(noexcept(std::declval<const view&>().num_bins()));
        static_assert(noexcept(std::declval<const view&>().size()));

        // The whole view is usable in a constant expression: construct it over
        // a compile-time buffer and read every accessor back.
        constexpr bool reads = [] {
            const float a[8] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};
            const view  s(a, 8);
            return s.size() == 8 && s.num_bins() == 5 && s.dc() == 1.0f && s.nyquist() == 2.0f && s.re(1) == 3.0f
                   && s.im(1) == 4.0f && s.power(0) == 1.0f && s.power(2) == 5.0f * 5.0f + 6.0f * 6.0f
                   && s.power(4) == 4.0f && s.bin_engineering(3).real() == 7.0f && s.bin_engineering(3).imag() == -8.0f;
        }();
        static_assert(reads);
        SUCCEED();
    }

} // namespace
