/// @file math.h
/// @brief The periodic Hann window and the decibel helpers, as consumer API.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The scalar formulas the tap::dsp primitives build on, promoted from
// tap/dsp/detail/math.h (Stage 6 of docs/audit-fft-and-code-smells.md) so the
// consuming libraries call the same functions instead of re-typing the
// expressions. tap::dsp::detail re-exports these three by using-declaration:
// there is one definition of each, and the names in the two namespaces are
// the same function.
//
// Contract (every point below is a contract point; changing one is a breaking
// numeric change for every consumer, and for pvoc.h, log_mel.h and the
// analysis instruments inside this repo):
//
//   periodic_hann(i, n) = 0.5 - 0.5 * std::cos(2.0 * pi * double(i) / double(n))
//   power_db(r)         = 10.0 * std::log10(r)
//   amplitude_db(r)     = 20.0 * std::log10(r)
//
// each evaluated in double, left to right, exactly as written: a call site
// that replaces the literal expression above with the helper moves no bit of
// its own arithmetic (tests/test_math.cpp pins that over a sweep; the
// same-host A/B fingerprint tool, tools/fingerprint, pins it for pvoc, log_mel
// and psola). One qualifier: under cross-statement contraction
// (-ffp-contract=fast, GCC's GNU-mode default) on FMA targets, any change to
// a translation unit's inlining decisions, this one included, can move OTHER
// code in that TU; read an A/B at those flags as a codegen change, and prove
// arithmetic identity at -ffp-contract=off (see the tools/fingerprint header
// and fft.h's D9 paragraph).
//
// Pi: this header exports no pi constant. Use std::numbers::pi (<numbers>),
// which is the double these helpers use; tap::dsp::detail::k_pi is the same
// double under another name and remains an implementation detail.

#pragma once

#include <cmath>
#include <cstddef>
#include <numbers>

namespace tap::dsp {

    /// The periodic Hann window, w[i] = 0.5 - 0.5 cos(2 pi i / n), i in [0, n):
    /// the n-point window whose n-periodic extension overlap-adds to a constant
    /// (the DFT-even form). It is NOT the symmetric window with denominator
    /// n - 1; a call site that uses the symmetric form must keep its own
    /// expression, since substituting this one changes every sample.
    ///
    /// Numerics (a contract point): evaluated in double as
    ///   0.5 - 0.5 * std::cos(((2.0 * std::numbers::pi) * i) / n)
    /// left to right, with i and n converted to double exactly (both below
    /// 2^53). w[0] is exactly 0.0; for even n, w[n/2] is exactly 1.0. The
    /// result is a double; a float profile narrows it once, at its own call
    /// site (as pvoc.h does), never by evaluating this in float.
    /// @pre n >= 1 and i < n.
    inline double periodic_hann(std::size_t i, std::size_t n) noexcept {
        return 0.5 - 0.5 * std::cos(2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(n));
    }

    /// Power ratio in decibels: 10 log10(ratio).
    ///
    /// Numerics (a contract point): evaluated as 10.0 * std::log10(ratio) in
    /// double. ratio == 0 gives -inf and ratio < 0 gives NaN, as std::log10
    /// does; a caller that needs a floor applies it to the argument
    /// (power_db(std::max(p, floor))), which keeps the association here.
    inline double power_db(double ratio) noexcept {
        return 10.0 * std::log10(ratio);
    }

    /// Amplitude ratio in decibels: 20 log10(ratio).
    ///
    /// Numerics (a contract point): evaluated as 20.0 * std::log10(ratio) in
    /// double, with the same edge behaviour as power_db.
    inline double amplitude_db(double ratio) noexcept {
        return 20.0 * std::log10(ratio);
    }

} // namespace tap::dsp
