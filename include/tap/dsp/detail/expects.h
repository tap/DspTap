/// @file expects.h
/// @brief TAP_EXPECTS — the house precondition check (STYLE.md §4), in its minimal form.
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

#pragma once

#include <cassert>

#if !defined(TAP_EXPECTS)
#define TAP_EXPECTS(cond) assert(cond)
#endif
