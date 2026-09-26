// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// COMPILE-FAIL FIXTURE, not a test executable: this translation unit must NOT
// compile. tests/CMakeLists.txt builds it as an excluded object library and
// the ctest `fft_compile_fail.LegacySpellingIsRejectedWithTheD4Message`
// (tests/compile_fail/expect_compile_failure.cmake) requires that the build
// fails, that the D4 message ("pre-Stage-4 spelling") is in the output, and
// that it is the ONLY error. fft.h promises that basic_real_fft<float |
// double, scaling::fixed> is rejected by name, first and alone (the D4
// expiry, tap/DspTap#40); without the named static_assert the concept assert
// would still reject it, with another message, so a plain WILL_FAIL would pass
// vacuously. The cv-qualified form is the same spelling and must take the same
// route.
#include "tap/dsp/fft.h"

#ifndef TAP_DSP_COMPILE_FAIL_CASE
#define TAP_DSP_COMPILE_FAIL_CASE 0
#endif

#if TAP_DSP_COMPILE_FAIL_CASE == 0
using legacy = tap::dsp::basic_real_fft<float, tap::dsp::scaling::fixed>;
#elif TAP_DSP_COMPILE_FAIL_CASE == 1
using legacy = tap::dsp::basic_real_fft<double, tap::dsp::scaling::fixed>;
#else
using legacy = tap::dsp::basic_real_fft<float, const tap::dsp::scaling::fixed>;
#endif

// Instantiates the class: the one thing the spelling may no longer do.
std::size_t tap_dsp_compile_fail_probe() {
    legacy fft(64);
    return fft.size();
}
