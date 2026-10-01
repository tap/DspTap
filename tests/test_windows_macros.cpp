// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Macro hygiene: every public header parses under the object-like macros the
// Windows SDK defines over ordinary names, so a consumer that includes
// <windows.h> before DspTap (the Max SDK does, under min-api) still compiles.
// The srdif engine once named locals `small`, which rpcndr.h defines as
// `char`; MSVC then rejected srdif.h in MuTap-Max's externals while every
// DspTap and MuTap leg, none of which includes <windows.h>, stayed green.
//
// The macros are defined here by hand, with the SDK's replacement lists
// (hyper's __int64 spelled portably), so the check runs on every host and on
// the QEMU legs rather than on the Windows leg alone; <windows.h> itself is
// not included. The standard library and gtest are included first, as a
// Windows consumer's would already be parsed: only DspTap's own text must
// survive the macros. min and max are out of scope: a consumer that includes
// <windows.h> defines NOMINMAX, as the SDK documents. Runs on every host and
// on all four QEMU legs (not excluded by MAIN_FILTER).

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <numbers>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

// The SDK's spellings, which the house macro-naming rule cannot apply to.
// NOLINTBEGIN(readability-identifier-naming)
// rpcndr.h
#ifndef small
#define small char
#endif
#ifndef hyper
#define hyper long long
#endif
// minwindef.h
#ifndef near
#define near
#endif
#ifndef far
#define far
#endif
// combaseapi.h
#ifndef interface
#define interface struct
#endif
// NOLINTEND(readability-identifier-naming)

#include "tap/dsp/analysis/multitone_analysis.h"
#include "tap/dsp/analysis/sine_analysis.h"
#include "tap/dsp/decimate.h"
#include "tap/dsp/fft.h"
#include "tap/dsp/fft/spectrum.h"
#include "tap/dsp/fir_kernels.h"
#include "tap/dsp/kaiser.h"
#include "tap/dsp/log_mel.h"
#include "tap/dsp/math.h"
#include "tap/dsp/nn.h"
#include "tap/dsp/psola.h"
#include "tap/dsp/pvoc.h"
#include "tap/dsp/quantize.h"
#include "tap/dsp/sample_traits.h"
#include "tap/dsp/yin.h"

namespace {

    // The engine's templates are instantiated as well as parsed: a round trip
    // through the default float and double engines at a size that runs the
    // fused small-block pass (N = 128, M = 64).
    template <typename Sample>
    Sample round_trip_error() {
        constexpr std::size_t            n = 128;
        tap::dsp::basic_real_fft<Sample> fft(n);
        std::vector<Sample>              a(n);
        for (std::size_t i = 0; i < n; ++i) {
            a[i] = static_cast<Sample>(std::sin(0.37 * static_cast<double>(i)));
        }
        const std::vector<Sample> x = a;
        fft.forward_inplace(a.data());
        fft.inverse_inplace(a.data());
        Sample worst = 0;
        for (std::size_t i = 0; i < n; ++i) {
            worst = std::max(worst, std::abs(a[i] * static_cast<Sample>(2.0 / n) - x[i]));
        }
        return worst;
    }

} // namespace

TEST(WindowsMacros, PublicHeadersParseUnderTheSdkMacros) {
    EXPECT_LT(round_trip_error<float>(), 1e-5F);
    EXPECT_LT(round_trip_error<double>(), 1e-13);
}
