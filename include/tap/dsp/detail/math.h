/// @file math.h
/// @brief Shared scalar math for the primitives: pi, plus the detail:: names of the public helpers.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// One home for the three formulas that were re-typed across the headers and
// tests (docs/audit-fft-and-code-smells.md, Part 2 "Duplication with no shared
// home"; Stage 6). Each helper is the exact expression its call sites already
// evaluated, in the same association order, so moving a call site here moves
// no output bit of that call site's own arithmetic: a helper call is one
// function boundary around the identical sequence of IEEE operations (the
// same-host A/B fingerprint tool, tools/fingerprint, is the gate that says so
// for pvoc, log_mel and psola). One qualifier: under GCC's cross-statement
// contraction (-ffp-contract=fast) on FMA targets, any change in a TU's
// inlining decisions, this one included, can move OTHER code in the same TU
// (see the tools/fingerprint header and fft.h's D9). Changing an association
// in tap/dsp/math.h is a numeric change for every consumer.
//
// Implementation detail of the tap::dsp headers; not a consumer-facing API.
// The periodic Hann window and the two dB helpers are defined in the public
// tap/dsp/math.h (their contract lives there) and re-exported here; k_pi stays
// here, since consumers take pi from std::numbers::pi directly.

#pragma once

#include <numbers>

#include "tap/dsp/math.h"

namespace tap::dsp::detail {

    /// pi as a double: std::numbers::pi. The literal 3.14159265358979323846
    /// that pvoc.h and psola.h carried rounds to the same double; the
    /// static_assert is the proof that replacing it changes no bit.
    inline constexpr double k_pi = std::numbers::pi;
    static_assert(k_pi == 3.14159265358979323846, "the retired literal and std::numbers::pi are one double");

    /// The periodic Hann window and the dB helpers are public since they were
    /// promoted to tap/dsp/math.h for consumers; these using-declarations keep
    /// the detail:: spellings the headers and tests use naming the same three
    /// functions (one definition each, so the two spellings cannot drift).
    using tap::dsp::amplitude_db;
    using tap::dsp::periodic_hann;
    using tap::dsp::power_db;

} // namespace tap::dsp::detail
