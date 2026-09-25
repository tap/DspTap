// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Contract battery for the public scalar helpers in tap/dsp/math.h: each one
// is, bit for bit, the literal expression its docstring states, so a consumer
// that replaces its own copy of that expression with the helper moves no bit.
// The inputs pass through a volatile first so neither side can be folded at
// compile time (GCC folds cos/log10 with MPFR, correctly rounded, which need
// not equal the run-time libm): both sides are evaluated by the same libm on
// the same operands, which is the claim a consumer's migration relies on.
// Runs on every host and on all four QEMU legs (not excluded by MAIN_FILTER).

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <type_traits>

#include <gtest/gtest.h>

#include "tap/dsp/detail/math.h"
#include "tap/dsp/math.h"

namespace {

    // Laundered through a volatile so the compiler cannot see the value.
    template <typename T>
    T opaque(T v) noexcept {
        volatile T sink = v;
        return sink;
    }

    std::uint64_t bits(double v) noexcept {
        return std::bit_cast<std::uint64_t>(v);
    }

    // The window lengths a consumer builds: every n up to 64 (all i), then the
    // frame sizes in use across the Tap libraries (odd and non-power-of-two
    // among them) and the powers of two to 4096.
    constexpr std::size_t k_lengths[] = {65,  96,  100, 127, 128,  160,  255,  256,  320,  400,
                                         480, 512, 640, 960, 1000, 1024, 1536, 2048, 4095, 4096};

    // One window length, every index: the helper against the literal.
    int count_hann_mismatches(std::size_t n) {
        int mismatches = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t ii      = opaque(i);
            const std::size_t nn      = opaque(n);
            const double      literal = 0.5 - 0.5 * std::cos(2.0 * std::numbers::pi * double(ii) / double(nn));
            if (bits(tap::dsp::periodic_hann(ii, nn)) != bits(literal)) {
                ++mismatches;
                ADD_FAILURE() << "periodic_hann(" << i << ", " << n << ") != literal";
            }
        }
        return mismatches;
    }

    // Ratios spanning what the batteries and instruments feed the dB helpers:
    // decades from 1e-300 to 1e300, a subnormal, exact powers of ten, and a
    // fixed-seed xorshift32 sweep of mantissas over [1e-12, 1e12).
    template <typename F>
    void for_each_ratio(F&& f) {
        for (int e = -300; e <= 300; e += 7) {
            f(std::pow(10.0, e));
        }
        f(std::numeric_limits<double>::denorm_min());
        f(std::numeric_limits<double>::min());
        f(std::numeric_limits<double>::max());
        for (const double r : {1.0, 2.0, 0.5, 10.0, 100.0, 0.1, 0.01, 3.0, 1e-6, 1e6}) {
            f(r);
        }
        std::uint32_t s = 0x9e3779b9U;
        for (int k = 0; k < 4096; ++k) {
            s ^= s << 13;
            s ^= s >> 17;
            s ^= s << 5;
            const double u = static_cast<double>(s) / 4294967296.0; // [0, 1)
            f(std::pow(10.0, -12.0 + 24.0 * u));
        }
    }

} // namespace

TEST(math, PeriodicHannIsTheLiteralExpressionBitForBit) {
    int mismatches = 0;
    for (std::size_t n = 1; n <= 64; ++n) {
        mismatches += count_hann_mismatches(n);
    }
    for (const std::size_t n : k_lengths) {
        mismatches += count_hann_mismatches(n);
    }
    EXPECT_EQ(mismatches, 0);
}

TEST(math, PeriodicHannIsZeroAtTheOriginAndOneAtTheCentre) {
    for (std::size_t n = 1; n <= 64; ++n) {
        EXPECT_EQ(bits(tap::dsp::periodic_hann(opaque(std::size_t{0}), opaque(n))), bits(0.0)) << "n = " << n;
    }
    for (std::size_t n = 2; n <= 4096; n += 2) {
        EXPECT_EQ(tap::dsp::periodic_hann(opaque(n / 2), opaque(n)), 1.0) << "n = " << n;
    }
}

TEST(math, PowerDbIsTheLiteralExpressionBitForBit) {
    int mismatches = 0;
    for_each_ratio([&](double r) {
        const double x = opaque(r);
        if (bits(tap::dsp::power_db(x)) != bits(10.0 * std::log10(x))) {
            ++mismatches;
            ADD_FAILURE() << "power_db(" << r << ") != 10.0 * std::log10";
        }
    });
    EXPECT_EQ(mismatches, 0);
}

TEST(math, AmplitudeDbIsTheLiteralExpressionBitForBit) {
    int mismatches = 0;
    for_each_ratio([&](double r) {
        const double x = opaque(r);
        if (bits(tap::dsp::amplitude_db(x)) != bits(20.0 * std::log10(x))) {
            ++mismatches;
            ADD_FAILURE() << "amplitude_db(" << r << ") != 20.0 * std::log10";
        }
    });
    EXPECT_EQ(mismatches, 0);
}

TEST(math, DbEdgesFollowLog10) {
    EXPECT_EQ(tap::dsp::power_db(opaque(1.0)), 0.0);
    EXPECT_EQ(tap::dsp::amplitude_db(opaque(1.0)), 0.0);
    EXPECT_EQ(tap::dsp::power_db(opaque(0.0)), -std::numeric_limits<double>::infinity());
    EXPECT_EQ(tap::dsp::amplitude_db(opaque(0.0)), -std::numeric_limits<double>::infinity());
    EXPECT_TRUE(std::isnan(tap::dsp::power_db(opaque(-1.0))));
    EXPECT_TRUE(std::isnan(tap::dsp::amplitude_db(opaque(-1.0))));
}

// The detail:: spellings the primitives use are using-declarations of the
// public functions: one definition each, so they cannot drift apart. Two
// non-type template arguments are the same type only when they name the same
// function.
namespace {
    template <auto F>
    struct function_identity {};
} // namespace
static_assert(
    std::is_same_v<function_identity<&tap::dsp::detail::periodic_hann>, function_identity<&tap::dsp::periodic_hann>>);
static_assert(std::is_same_v<function_identity<&tap::dsp::detail::power_db>, function_identity<&tap::dsp::power_db>>);
static_assert(
    std::is_same_v<function_identity<&tap::dsp::detail::amplitude_db>, function_identity<&tap::dsp::amplitude_db>>);
static_assert(tap::dsp::detail::k_pi == std::numbers::pi);

TEST(math, DetailNamesAgreeWithThePublicNames) {
    for (const std::size_t n : k_lengths) {
        for (std::size_t i = 0; i < n; i += 7) {
            EXPECT_EQ(bits(tap::dsp::detail::periodic_hann(opaque(i), opaque(n))),
                      bits(tap::dsp::periodic_hann(opaque(i), opaque(n))));
        }
    }
    for_each_ratio([](double r) {
        EXPECT_EQ(bits(tap::dsp::detail::power_db(opaque(r))), bits(tap::dsp::power_db(opaque(r))));
        EXPECT_EQ(bits(tap::dsp::detail::amplitude_db(opaque(r))), bits(tap::dsp::amplitude_db(opaque(r))));
    });
}
