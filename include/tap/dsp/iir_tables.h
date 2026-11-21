/// @file iir_tables.h
/// @brief The committed recursive designs: near-linear half-bands, four-path stages, Chebyshev II gains, type II FIRs.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The designs iir_design.h cannot compute at compile time, committed as
// tables (SampleRateTap's docs/IIR_PLAN.md, I4): each comes from an iterative
// optimization (a linear program and SLSQP, a minimax over random starts,
// Parks-McClellan) or a sweep too long for MSVC's constant evaluator (the
// L-inf gains), so tools/iir/gen_iir_tables.py generates them, in a venv
// with the versions tools/iir/requirements.txt pins, into
// detail/iir_tables_data.h, and the tests pin what each row promises
// (tests/test_iir_design.cpp). Regenerate only when a method or a spec
// point changes, and commit the data header with the change.
//
// Every row is keyed by its spec point, as numbers (iir_design.h's file
// header): stop_db, the passband edge pass_num / pass_den of the chain's
// lowest rate, and a key, the design divisor for the half-bands and the
// four-path stage, the Chebyshev order for the guard band filter's rows.
// iir_design.h looks a row up by spec at compile time
// (k_near_linear_halfband<...>, k_fourpath<...>, k_guard_band_cheby2<...>,
// k_guard_band_fir<...>); a missing row is a static_assert.

#pragma once

#include <array>
#include <cstddef>
#include <span>

namespace tap::dsp {

    /// The most sections a Chebyshev II design holds: orders up to 15.
    inline constexpr std::size_t k_max_cheby2_sections = 8;

    /// A near-linear two-path half-band H(z) = 1/2 [z^-2K + z^-1 A1(z^2)]:
    /// path 0 the pure delay of K path (low-rate) samples, path 1 the allpass
    /// of order N in w = z^2, as sections: (c1, c2) per second-order section
    /// (c2 + c1 w^-1 + w^-2) / (1 + c1 w^-1 + c2 w^-2), ordered by increasing
    /// pole radius, then c1 alone for one first-order section (c1 + w^-1) /
    /// (1 + c1 w^-1) when N is odd.
    struct near_linear_halfband_design {
        int                     stop_db  = 0;
        int                     pass_num = 0;
        int                     pass_den = 1;
        int                     divisor  = 1;
        int                     order    = 0; ///< N: multiplies per decimated output
        int                     delay    = 0; ///< K, path samples
        std::span<const double> sections;
        double                  stopband_db     = 0.0; ///< attained, on 16384 points of [1/2 - f_pass, 1/2]
        double                  delay_deviation = 0.0; ///< max |tau(f) - 2K| over the passband, high-rate samples
        double                  pole_radius     = 0.0; ///< the largest, in path samples

        constexpr int key() const noexcept { return divisor; }
    };

    /// A four-path stage H(z) = 1/4 sum_k z^-k A_k(z^4), each A_k a cascade of
    /// first-order allpass sections (a + z^-4) / (1 + a z^-4).
    struct fourpath_design {
        int                                    stop_db  = 0;
        int                                    pass_num = 0;
        int                                    pass_den = 1;
        int                                    divisor  = 1; ///< the stage's lower rate over the chain's lowest
        std::array<std::span<const double>, 4> paths{};
        double stopband_db = 0.0; ///< attained on the bands that fold onto the passband (16384 points)

        constexpr int key() const noexcept { return divisor; }
    };

    /// The L-inf gains of a Chebyshev II designed at unity DC gain per
    /// section: gains[i] scales section i's numerator, i < count - 1 (the
    /// last section takes 1 / their product).
    struct cheby2_gains {
        int                                           stop_db  = 0;
        int                                           pass_num = 0;
        int                                           pass_den = 1;
        int                                           order    = 0;
        std::array<double, k_max_cheby2_sections - 1> gains{};

        constexpr int key() const noexcept { return order; }
    };

    /// A linear-phase (type II, even-length, symmetric) guard band FIR, flat
    /// within 0.01 dB peak to peak to its passband and stop_db down from the
    /// order-`order` Chebyshev II's stop edge at the same spec point.
    struct guard_band_fir_design {
        int                     stop_db   = 0;
        int                     pass_num  = 0;
        int                     pass_den  = 1;
        int                     order     = 0;   ///< the Chebyshev II whose stop edge it matches
        double                  stop_edge = 0.0; ///< fraction of its rate
        std::span<const double> taps;
        double                  ripple_db   = 0.0; ///< passband, peak to peak, attained
        double                  stopband_db = 0.0; ///< attained, from stop_edge

        constexpr int key() const noexcept { return order; }
    };

} // namespace tap::dsp

#include "tap/dsp/detail/iir_tables_data.h"
