/// @file math.h
/// @brief Shared scalar math for the primitives: pi, the periodic Hann window, decibels.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// One home for the three formulas that were re-typed across the headers and
// tests (docs/audit-fft-and-code-smells.md, Part 2 "Duplication with no shared
// home"; Stage 6). Each helper is the exact expression its call sites already
// evaluated, in the same association order, so moving a call site here moves
// no output bit: a helper call is one function boundary around the identical
// sequence of IEEE operations (the same-host A/B fingerprint tool,
// tools/fingerprint, is the gate that says so for pvoc, log_mel and psola).
// Changing an association below is a numeric change for every consumer.
//
// Implementation detail of the tap::dsp headers; not a consumer-facing API.

#pragma once

#include <cmath>
#include <cstddef>
#include <numbers>

namespace tap::dsp::detail {

    /// pi as a double: std::numbers::pi. The literal 3.14159265358979323846
    /// that pvoc.h and psola.h carried rounds to the same double; the
    /// static_assert is the proof that replacing it changes no bit.
    inline constexpr double k_pi = std::numbers::pi;
    static_assert(k_pi == 3.14159265358979323846, "the retired literal and std::numbers::pi are one double");

    /// The periodic Hann window, w[i] = 0.5 - 0.5 cos(2 pi i / n), i in [0, n):
    /// the n-point window whose n-periodic extension overlap-adds to a constant
    /// (the DFT-even form; not the symmetric n - 1 denominator).
    ///
    /// Numerics (a contract point): evaluated in double as
    ///   0.5 - 0.5 * std::cos(((2.0 * pi) * i) / n)
    /// left to right, with i and n converted to double exactly (both below
    /// 2^53). pvoc.h and log_mel.h build their windows from this expression and
    /// their outputs are bit-pinned to it by the fingerprint tool.
    /// @pre n >= 1 and i < n.
    inline double periodic_hann(std::size_t i, std::size_t n) noexcept {
        return 0.5 - 0.5 * std::cos(2.0 * k_pi * static_cast<double>(i) / static_cast<double>(n));
    }

    /// Power ratio in decibels: 10 log10(ratio). Evaluated as
    /// 10.0 * std::log10(ratio); ratio <= 0 gives -inf / NaN as std::log10 does
    /// (callers that need a floor apply it to the argument).
    inline double power_db(double ratio) noexcept {
        return 10.0 * std::log10(ratio);
    }

    /// Amplitude ratio in decibels: 20 log10(ratio). Evaluated as
    /// 20.0 * std::log10(ratio), with the same edge behaviour as power_db.
    inline double amplitude_db(double ratio) noexcept {
        return 20.0 * std::log10(ratio);
    }

} // namespace tap::dsp::detail
