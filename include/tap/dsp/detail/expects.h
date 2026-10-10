/// @file expects.h
/// @brief The two contract idioms: TAP_EXPECTS, the house precondition check (STYLE.md §4) in
/// its minimal form, and tap::dsp::raise, the construction-time rejection that throws where
/// the build has exceptions and terminates where it does not.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// TAP_EXPECTS(cond) states a precondition the way STYLE.md §4 prescribes for
// the Tap libraries: an assertion in a debug build (NDEBUG not defined) that
// evaluates to nothing in a release build, and never throws — some Tap
// consumers build with -fno-exceptions, and a precondition is the caller's
// bug, not a recoverable condition. It is `assert` under a name that says
// "contract", so the sites that state a documented @pre are greppable and a
// future repo-wide precondition policy (its own plan; audit Stage 6 lands
// this macro for fft.h's power-of-two and size-range preconditions only, at
// Stage 4) can change what every site does in one place.
//
// What it does NOT do, stated so nobody relies on it: it does not check in a
// release build. A violated precondition in a release build is undefined
// behaviour exactly as it was under the plain assert it replaces; the
// release-mode tool is the constexpr predicate the class offers beside the
// precondition (basic_real_fft<>::supports_size), which a caller checks
// before constructing. Every CI battery runs Release / MinSizeRel, so a test
// of the assertion itself runs only in a local Debug configure.
//
// A consumer may define TAP_EXPECTS before including any tap/dsp header to
// substitute its own policy (a release-mode terminate, a logging hook). The
// definition must be the same in every translation unit of an image: the
// macro expands inside inline and template functions, and two definitions
// in one image are an ODR violation. Across images the same shape recurs
// as the cross-image hazard fft.h describes for layouts (audit F4), here
// code-only: two images that define TAP_EXPECTS differently and both
// instantiate the same inline function export two weak definitions under
// one name, and a loader that coalesces them makes one image run the
// other's precondition policy. Same remedy: keep the definition uniform
// across the images of one process, or hide visibility.
//
// tap::dsp::raise<E>(args...) is for the other kind of failure: not a
// precondition the caller is bound to meet, but a construction-time
// rejection of a value the caller could not have checked without the
// library's own arithmetic — a configuration that fails the design rules,
// an allocation the platform refused. Where the translation unit is
// compiled with C++ exceptions (TAP_DSP_HAS_EXCEPTIONS is 1: __cpp_exceptions
// on GCC and Clang, _CPPUNWIND under MSVC's /EH) it is `throw E(args...)`,
// the contract every hosted consumer and test battery pins
// (std::invalid_argument on a bad configuration, std::bad_alloc on a refused
// allocation). Compiled with exceptions off (-fno-exceptions, the default of
// the Pico SDK, Zephyr and most Cortex-M firmware trees) a `throw` does not
// compile and no handler for E exists anywhere in the image, so the call is
// a contract violation and terminates — std::terminate(), through the
// handler the firmware installed, if any — exactly as the standard
// library's own allocation failure already does there. Every header that
// rejects through raise therefore compiles either way, and the release-mode
// tool is the same as for preconditions: the predicate the class offers
// beside the rejection (a validate() that returns the reason, a
// supports_size), checked before constructing, so that such a firmware never
// reaches the call. The ODR hazard above applies unchanged: raise expands
// inside inline and template functions, so every translation unit of an
// image must agree on TAP_DSP_HAS_EXCEPTIONS (one exception setting per
// image).

#pragma once

#include <cassert>
#include <exception>
#include <utility>

#if !defined(TAP_EXPECTS)
#define TAP_EXPECTS(cond) assert(cond)
#endif

#if defined(__cpp_exceptions) || defined(_CPPUNWIND)
#define TAP_DSP_HAS_EXCEPTIONS 1
#else
#define TAP_DSP_HAS_EXCEPTIONS 0
#endif

namespace tap::dsp {

    /// Rejects at construction time: throws E(args...) where the build has
    /// exceptions, std::terminate()s where it does not (the file comment).
    template <typename E, typename... Args>
    [[noreturn]] inline void raise(Args&&... args) {
#if TAP_DSP_HAS_EXCEPTIONS
        throw E(std::forward<Args>(args)...);
#else
        (static_cast<void>(args), ...);
        std::terminate();
#endif
    }

} // namespace tap::dsp
