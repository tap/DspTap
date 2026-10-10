// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Pins the two contract idioms of tap/dsp/detail/expects.h as the battery
// builds them. Every test target compiles with exceptions (gtest's
// EXPECT_THROW needs them), so TAP_DSP_HAS_EXCEPTIONS is 1 here and
// raise<E>(args...) must throw exactly E, constructed from the caller's
// arguments. The exceptions-off side — raise terminates, every public header
// compiles — is the compile check beside this file in tests/CMakeLists.txt
// (expects.PublicHeadersCompileWithoutExceptions): a death test of the
// terminate would need an exceptions-off gtest, which the battery does not
// build. TAP_EXPECTS is pinned for the one property a release build can
// observe: under NDEBUG its argument is not evaluated (the CI batteries run
// Release / MinSizeRel; a local Debug configure sees the assert instead).

#include <new>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "tap/dsp/detail/expects.h"

namespace {

    TEST(Expects, TheBatteryBuildsWithExceptions) {
        static_assert(TAP_DSP_HAS_EXCEPTIONS == 1, "the test battery needs exceptions (EXPECT_THROW)");
        EXPECT_EQ(TAP_DSP_HAS_EXCEPTIONS, 1);
    }

    TEST(Expects, RaiseThrowsTheExceptionItNamesWithTheCallersText) {
        EXPECT_THROW(tap::dsp::raise<std::invalid_argument>("bad configuration"), std::invalid_argument);
        try {
            tap::dsp::raise<std::invalid_argument>("tap::dsp::raise: the caller's text");
            ADD_FAILURE() << "raise returned";
        }
        catch (const std::invalid_argument& e) {
            EXPECT_STREQ(e.what(), "tap::dsp::raise: the caller's text");
        }
        // A std::string argument is forwarded as the exception's constructor
        // takes it, so a composed reason needs no c_str().
        const std::string reason = "composed";
        EXPECT_THROW(tap::dsp::raise<std::runtime_error>(reason), std::runtime_error);
    }

    TEST(Expects, RaiseWithoutArgumentsThrowsADefaultConstructedException) {
        // accelerate.h's refused vDSP table: bad_alloc has no message.
        EXPECT_THROW(tap::dsp::raise<std::bad_alloc>(), std::bad_alloc);
    }

    TEST(Expects, RaiseIsNoreturn) {
        // Compiles only because raise is [[noreturn]]: the function has no
        // return statement on this path and -Wreturn-type is an error in the
        // warnings target.
        const auto f = [](bool reject) -> int {
            if (reject) {
                tap::dsp::raise<std::invalid_argument>("rejected");
            }
            return 7;
        };
        EXPECT_EQ(f(false), 7);
        EXPECT_THROW(f(true), std::invalid_argument);
    }

    TEST(Expects, ExpectsDoesNotEvaluateItsArgumentInRelease) {
        int        evaluations = 0;
        const auto cond        = [&evaluations]() {
            ++evaluations;
            return true;
        };
        TAP_EXPECTS(cond());
        static_cast<void>(cond); // referenced only inside the macro, which NDEBUG drops
#if defined(NDEBUG)
        EXPECT_EQ(evaluations, 0) << "a release build's precondition costs nothing";
#else
        EXPECT_EQ(evaluations, 1) << "a debug build's precondition is the assert";
#endif
    }

} // namespace
