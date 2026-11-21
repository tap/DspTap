// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Pins subnormals.h: which sample types a control bit governs on this target,
// that subnormals_flushed<S>() reads the bit a host sets and clears, that it
// never writes the control word, and that it is false wherever no bit exists.

#include <cstdint>
#include <limits>

#include <gtest/gtest.h>

#include "support/flush_to_zero.h"
#include "tap/dsp/subnormals.h"

namespace {

    using tap::dsp::has_subnormal_flush_control;
    using tap::dsp::subnormals_flushed;
    using tap_test::fp_control_bits;
    using tap_test::scoped_flush_to_zero;

    TEST(Subnormals, TheTraitIsTheTargetTable) {
        static_assert(!has_subnormal_flush_control<std::int16_t>());
        static_assert(!has_subnormal_flush_control<std::int32_t>());
#if defined(__x86_64__) || defined(_M_X64) || defined(__aarch64__) || defined(_M_ARM64)
        static_assert(has_subnormal_flush_control<float>());
        static_assert(has_subnormal_flush_control<double>());
#elif defined(__ARM_FP)
        static_assert(has_subnormal_flush_control<float>() == ((__ARM_FP & 0x4) != 0));
        static_assert(has_subnormal_flush_control<double>() == ((__ARM_FP & 0x8) != 0));
#elif defined(__arm__) || defined(__hexagon__)
        static_assert(!has_subnormal_flush_control<float>());
        static_assert(!has_subnormal_flush_control<double>());
#endif
        SUCCEED();
    }

    template <typename S>
    void reads_the_bit() {
        if constexpr (has_subnormal_flush_control<S>()) {
            {
                scoped_flush_to_zero on(true);
                EXPECT_TRUE(subnormals_flushed<S>());
                // And the hardware agrees: a subnormal result is zero.
                volatile S tiny = std::numeric_limits<S>::min();
                volatile S half = tiny / S(2);
                EXPECT_EQ(half, S(0));
            }
            {
                scoped_flush_to_zero off(false);
                EXPECT_FALSE(subnormals_flushed<S>());
                volatile S tiny = std::numeric_limits<S>::min();
                volatile S half = tiny / S(2);
                EXPECT_NE(half, S(0));
            }
        }
        else {
            EXPECT_FALSE(subnormals_flushed<S>());
        }
    }

    TEST(Subnormals, ReadsTheBitAHostSetsFloat) {
        reads_the_bit<float>();
    }
    TEST(Subnormals, ReadsTheBitAHostSetsDouble) {
        reads_the_bit<double>();
    }

    TEST(Subnormals, NeverWritesTheControlWord) {
        for (const bool on : {true, false}) {
            scoped_flush_to_zero mode(on);
            const auto           before = fp_control_bits();
            static_cast<void>(subnormals_flushed<float>());
            static_cast<void>(subnormals_flushed<double>());
            EXPECT_EQ(fp_control_bits(), before);
        }
    }

    TEST(Subnormals, FixedPointIsNeverFlushed) {
        EXPECT_FALSE(subnormals_flushed<std::int16_t>());
        EXPECT_FALSE(subnormals_flushed<std::int32_t>());
    }

} // namespace
