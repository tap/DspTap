// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// tap/dsp/version.h: the header's version is the CMake project()'s (the
// build passes PROJECT_VERSION_* in as TAP_DSP_TEST_PROJECT_VERSION_*), and
// the encoding is the documented integer packing, so a consumer's
// static_assert floor compares the way it reads.

#include <gtest/gtest.h>

#include "tap/dsp/version.h"

namespace {

    TEST(Version, HeaderMatchesTheProjectVersion) {
        EXPECT_EQ(TAP_DSP_VERSION_MAJOR, TAP_DSP_TEST_PROJECT_VERSION_MAJOR);
        EXPECT_EQ(TAP_DSP_VERSION_MINOR, TAP_DSP_TEST_PROJECT_VERSION_MINOR);
        EXPECT_EQ(TAP_DSP_VERSION_PATCH, TAP_DSP_TEST_PROJECT_VERSION_PATCH);
    }

    TEST(Version, EncodingIsBitPackedAndOrdered) {
        static_assert(TAP_DSP_VERSION_ENCODE(0, 1, 0) == 0x000100);
        static_assert(TAP_DSP_VERSION_ENCODE(1, 2, 3) == 0x010203);
        static_assert(TAP_DSP_VERSION
                      == (TAP_DSP_VERSION_MAJOR << 16 | TAP_DSP_VERSION_MINOR << 8 | TAP_DSP_VERSION_PATCH));
        // A patch never outranks a minor, nor a minor a major.
        static_assert(TAP_DSP_VERSION_ENCODE(0, 1, 255) < TAP_DSP_VERSION_ENCODE(0, 2, 0));
        static_assert(TAP_DSP_VERSION_ENCODE(0, 255, 255) < TAP_DSP_VERSION_ENCODE(1, 0, 0));
#if TAP_DSP_VERSION < TAP_DSP_VERSION_ENCODE(0, 1, 0)
        FAIL() << "the encoded version must work in #if";
#endif
        SUCCEED();
    }

} // namespace
