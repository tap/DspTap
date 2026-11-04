// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Pins the two contract idioms of tap/dsp/detail/expects.h as the battery
// builds them, in both of the configurations it is built in.
//
// With exceptions (the default), TAP_DSP_HAS_EXCEPTIONS is 1 and
// raise<E>(args...) must throw exactly E, constructed from the caller's
// arguments. Built with TAP_DSP_NO_EXCEPTIONS_RTTI (the whole battery under
// -fno-exceptions -fno-rtti, CI's no-exceptions legs), it is 0 and raise must
// terminate: a death test where gtest has them (the hosted legs), and on bare
// metal, where a terminate ends the one-shot run, the property is that the
// build reached every translation unit (TheBatteryBuildsWithoutExceptionsOrRtti)
// and that raise is still [[noreturn]]. Every public header compiling under
// -fno-exceptions alone is the compile check beside this file in
// tests/CMakeLists.txt (expects.PublicHeadersCompileWithoutExceptions).
// TAP_EXPECTS is pinned for the one property a release build can observe:
// under NDEBUG its argument is not evaluated (the CI batteries run Release /
// MinSizeRel; a local Debug configure sees the assert instead).

#include <new>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "tap/dsp/detail/expects.h"

namespace {

#if TAP_DSP_HAS_EXCEPTIONS

    TEST(Expects, TheBatteryBuildsWithExceptions) {
        static_assert(TAP_DSP_HAS_EXCEPTIONS == 1, "this half of the file is the exceptions build");
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

#else

    // The option reached this translation unit: neither feature is on, so
    // every other file in the battery was compiled the same way.
    TEST(Expects, TheBatteryBuildsWithoutExceptionsOrRtti) {
        static_assert(TAP_DSP_HAS_EXCEPTIONS == 0, "this half of the file is the exceptions-off build");
#if defined(TAP_DSP_NO_EXCEPTIONS_RTTI)
#if defined(__cpp_rtti) || defined(__GXX_RTTI)
#error "TAP_DSP_NO_EXCEPTIONS_RTTI is set but RTTI is on"
#endif
#endif
        EXPECT_EQ(TAP_DSP_HAS_EXCEPTIONS, 0);
    }

#if GTEST_HAS_DEATH_TEST
    // Where the build has no exceptions a rejection is a contract violation:
    // raise terminates (std::terminate(), so the process dies on SIGABRT).
    TEST(ExpectsDeathTest, RaiseTerminatesWithoutExceptions) {
        EXPECT_DEATH(tap::dsp::raise<std::invalid_argument>("bad configuration"), "");
        const std::string reason = "composed";
        EXPECT_DEATH(tap::dsp::raise<std::runtime_error>(reason), "");
    }

    TEST(ExpectsDeathTest, RaiseWithoutArgumentsTerminatesWithoutExceptions) {
        EXPECT_DEATH(tap::dsp::raise<std::bad_alloc>(), "");
    }
#endif

#endif

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
#if TAP_DSP_HAS_EXCEPTIONS
        EXPECT_THROW(f(true), std::invalid_argument);
#endif // without exceptions, the rejecting path is ExpectsDeathTest's
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
