// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// THE SRDIF ENGINE'S OWN CONTRACT POINTS (include/tap/dsp/fft/srdif.h), beyond
// what every floating engine promises (test_fft.cpp, test_fft_oracle.cpp,
// test_fft_engine.cpp): the integer trigonometry its tables are built from.
//
//   - srdif_trig::first_octant evaluates cos, sin, 1 - cos and 1 - sin of
//     2 pi k / 2^L in 64-bit fixed point, and to_sample rounds that to the
//     nearest float / double. Against std::cos / std::sin of the double
//     angle (which carry their own <= 1 ulp and the angle's rounding) the
//     float values agree to within one float ulp and the double values to
//     within the pinned ulp count below; the precision statement proper
//     (every float entry correctly rounded, every double entry within
//     0.5 + 2^-7.3 ulp) is a quad-precision measurement recorded in the
//     header, not something a portable test can re-derive.
//   - to_sample rounds half to even on the 64-bit mantissa, and the
//     constants the kernel spells as literals (sqrt(1/2), cos and sin of
//     pi/8) are the generator's own values, so the leaves and the tables
//     agree bit for bit.
//
// The generator is integer arithmetic plus one exact conversion and one
// exact power-of-two scaling, so these values are the same on every host;
// the output fingerprints (test_fft_srdif_fingerprint.cpp) pin that
// end to end.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numbers>

#include <gtest/gtest.h>

#include "tap/dsp/fft/srdif.h"

namespace {

    namespace trig = tap::dsp::detail::srdif_trig;

    /// |a - b| in units of b's ulp (b != 0).
    template <typename Sample>
    double ulps(Sample a, Sample b) {
        const Sample ulp = std::nextafter(std::fabs(b), std::numeric_limits<Sample>::infinity()) - std::fabs(b);
        return static_cast<double>(std::fabs(a - b) / ulp);
    }

    // Measured 2026-09-26 on x86-64 glibc 2.39 (both libm dispatches) over
    // n = 2^3 … 2^14, every k of the first octant: cos and sin 0 float ulp
    // and 1 double ulp at most (libm's last bit and the double angle's
    // rounding); 1 - sin against the double expression 1 - std::sin, 2 ulp
    // (the expression's own rounding of sin, carried into a result up to
    // 3.4x smaller); 1 - cos against 2 sin^2(theta/2), 3 ulp (three
    // roundings of the expression; 4 ulp on the QEMU legs' newlib). Pinned
    // at 1 / 2 / 4 / 6 ulp: these are
    // sanity bounds against a second implementation, not the precision
    // statement (the file comment).
    constexpr double k_float_ulps         = 1.0;
    constexpr double k_double_ulps        = 2.0;
    constexpr double k_one_minus_sin_ulps = 4.0;
    constexpr double k_one_minus_cos_ulps = 6.0;

    TEST(fft_srdif_trig, GeneratorAgreesWithLibmToTheLastBits) {
        double worst_f   = 0.0;
        double worst_d   = 0.0;
        double worst_oms = 0.0;
        double worst_omc = 0.0;
        for (int lg = 3; lg <= 14; ++lg) {
            const std::uint64_t n = std::uint64_t{1} << lg;
            for (std::uint64_t k = 1; k <= n / 8; ++k) {
                const auto   v     = trig::first_octant(k, lg);
                const double angle = 2.0 * std::numbers::pi * static_cast<double>(k) / static_cast<double>(n);
                const double c     = std::cos(angle);
                const double s     = std::sin(angle);
                const double h     = std::sin(angle / 2.0);
                worst_d            = std::fmax(worst_d, ulps(trig::to_sample<double>(v.cos), c));
                worst_d            = std::fmax(worst_d, ulps(trig::to_sample<double>(v.sin), s));
                worst_f            = std::fmax(worst_f, ulps(trig::to_sample<float>(v.cos), static_cast<float>(c)));
                worst_f            = std::fmax(worst_f, ulps(trig::to_sample<float>(v.sin), static_cast<float>(s)));
                worst_oms          = std::fmax(worst_oms, ulps(trig::to_sample<double>(v.one_minus_sin), 1.0 - s));
                worst_omc          = std::fmax(worst_omc, ulps(trig::to_sample<double>(v.one_minus_cos), 2.0 * h * h));
            }
        }
        std::printf("[ measured ] srdif trig vs libm: cos/sin float %.2f ulp, double %.2f ulp; 1 - sin %.2f ulp, "
                    "1 - cos %.2f ulp\n",
                    worst_f, worst_d, worst_oms, worst_omc);
        EXPECT_LE(worst_f, k_float_ulps);
        EXPECT_LE(worst_d, k_double_ulps);
        EXPECT_LE(worst_oms, k_one_minus_sin_ulps);
        EXPECT_LE(worst_omc, k_one_minus_cos_ulps);
    }

    TEST(fft_srdif_trig, OctantEndpointsAreExact) {
        for (int lg = 3; lg <= 30; ++lg) {
            const auto zero = trig::first_octant(0, lg);
            EXPECT_EQ(trig::to_sample<double>(zero.cos), 1.0);
            EXPECT_EQ(trig::to_sample<double>(zero.sin), 0.0);
            EXPECT_EQ(trig::to_sample<double>(zero.one_minus_cos), 0.0);
            EXPECT_EQ(trig::to_sample<double>(zero.one_minus_sin), 1.0);
            // k = n/8: cos = sin = sqrt(1/2), rounded identically.
            const auto eighth = trig::first_octant(std::uint64_t{1} << (lg - 3), lg);
            EXPECT_EQ(trig::to_sample<double>(eighth.cos), 0x1.6a09e667f3bcdp-1) << "L=" << lg;
            EXPECT_EQ(trig::to_sample<double>(eighth.sin), 0x1.6a09e667f3bcdp-1) << "L=" << lg;
            EXPECT_EQ(trig::to_sample<float>(eighth.cos), 0x1.6a09e6p-1f) << "L=" << lg;
            EXPECT_EQ(trig::to_sample<float>(eighth.sin), 0x1.6a09e6p-1f) << "L=" << lg;
        }
    }

    // The literals the 16-point leaf and the fused pass spell out are the
    // generator's values at n = 16, k = 1.
    TEST(fft_srdif_trig, KernelLiteralsAreTheGeneratorsValues) {
        const auto v = trig::first_octant(1, 4);
        EXPECT_EQ(trig::to_sample<double>(v.cos), 0x1.d906bcf328d46p-1);
        EXPECT_EQ(trig::to_sample<double>(v.sin), 0x1.87de2a6aea963p-2);
        EXPECT_EQ(trig::to_sample<float>(v.cos), 0x1.d906bcp-1f);
        EXPECT_EQ(trig::to_sample<float>(v.sin), 0x1.87de2ap-2f);
    }

    TEST(fft_srdif_trig, ToSampleRoundsHalfToEven) {
        // m = 2^63 + t, value m * 2^-63 = 1 + t 2^-63. Float keeps 24 bits:
        // its ulp at 1 is 2^-23, i.e. t = 2^40 per ulp.
        const std::uint64_t one = std::uint64_t{1} << 63;
        const auto          f   = [&](std::uint64_t t) { return trig::to_sample<float>({one + t, -63}); };
        EXPECT_EQ(f(0), 1.0f);
        EXPECT_EQ(f((std::uint64_t{1} << 39) - 1), 1.0f);                  // below half: down
        EXPECT_EQ(f(std::uint64_t{1} << 39), 1.0f);                        // tie, even below: down
        EXPECT_EQ(f((std::uint64_t{1} << 39) + 1), 1.0f + 0x1p-23f);       // above half: up
        EXPECT_EQ(f((std::uint64_t{3} << 39)), 1.0f + 0x1p-22f);           // tie, odd below: up
        EXPECT_EQ(trig::to_sample<float>({~std::uint64_t{0}, -64}), 1.0f); // carries into the exponent
        // Double keeps 53 bits: ulp 2^-52 at 1, t = 2^11 per ulp.
        const auto d = [&](std::uint64_t t) { return trig::to_sample<double>({one + t, -63}); };
        EXPECT_EQ(d(std::uint64_t{1} << 10), 1.0);
        EXPECT_EQ(d((std::uint64_t{1} << 10) + 1), 1.0 + 0x1p-52);
        EXPECT_EQ(d(std::uint64_t{3} << 10), 1.0 + 0x1p-51);
        // Scaling is exact far from 1 (the smallest angle of the largest size).
        EXPECT_EQ(trig::to_sample<double>({one, -63 - 40}), 0x1p-40);
    }

} // namespace
