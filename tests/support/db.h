/// @file db.h
/// @brief The test battery's decibel helpers: the library's own, re-exported.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The tests convert to dB with the same two functions the instruments use
// (include/tap/dsp/detail/math.h), so a test and the header it checks cannot
// disagree about what "dB" means. power_db is 10 log10, amplitude_db 20 log10,
// both evaluated exactly as the inline expressions they replaced at Stage 6;
// a floor, where a test needs one, is applied to the argument at the call.

#pragma once

#include "tap/dsp/detail/math.h"

namespace tap::dsp::test {

    using tap::dsp::detail::amplitude_db;
    using tap::dsp::detail::power_db;

} // namespace tap::dsp::test
