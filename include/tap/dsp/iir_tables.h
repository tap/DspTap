/// @file iir_tables.h
/// @brief The committed recursive designs: near-linear half-bands, four-path stages, Chebyshev II gains, type II FIRs
/// and their Q15 rows.
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
// k_guard_band_fir<...>, k_guard_band_fir_q15<...>); a missing row is a
// static_assert.

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

    /// A type II guard band FIR designed for Q15 (DspTap 0.4.0): a double
    /// design whose taps, quantized to Q1.14 by
    /// quantize_symmetric_row_preserving_sum<std::int16_t> (quantize.h), keep
    /// the double row's spec: flat within 0.01 dB peak to peak to the passband
    /// and stop_db down from the same stop edge, the order-`order` Chebyshev
    /// II's. The double rows lose it in Q1.14: quantized the same way, the ones
    /// at 70 dB attain -63.5 to -74.8 dB from their stop edges, and thirteen of
    /// the fifteen miss the stopband or the ripple (tests/test_iir_design.cpp).
    /// So each Q15 row buys the margin quantization needs. It is the shortest
    /// even length (from 4) at which a Parks-McClellan design (McClellan, Parks
    /// and Rabiner 1973, through scipy.signal.remez), targeted at stop_db to
    /// stop_db + 30 dB in 0.01 dB steps, meets the spec in Q1.14 on
    /// 16384-point grids over the passband and over the stopband; at that
    /// length, the design whose Q1.14 taps attain the deepest stopband
    /// (tools/iir/gen_iir_tables.py). In remez only the ratio of the band
    /// weights matters, so a higher target and more stopband weight are one
    /// lever, and the target spans both. Rows at the three 70 dB spec points
    /// and the double table's orders 3 to 11, the stopband and ripple those of
    /// the Q1.14 taps:
    ///
    ///   | A, p          | order | taps (double row) | stopband  | ripple     | designed for |
    ///   |---------------|-------|-------------------|-----------|------------|--------------|
    ///   | 70 dB, 1/3    | 3     | 20 (20)           | -75.52 dB | 0.00965 dB | 73.65 dB     |
    ///   | 70 dB, 1/3    | 5     | 32 (30)           | -89.23 dB | 0.00898 dB | 89.27 dB     |
    ///   | 70 dB, 1/3    | 7     | 40 (38)           | -70.40 dB | 0.00984 dB | 70.79 dB     |
    ///   | 70 dB, 1/3    | 9     | 54 (52)           | -72.82 dB | 0.00939 dB | 74.40 dB     |
    ///   | 70 dB, 1/3    | 11    | 72 (66)           | -73.12 dB | 0.00982 dB | 74.27 dB     |
    ///   | 70 dB, 3/8    | 3     | 26 (26)           | -76.45 dB | 0.00970 dB | 75.46 dB     |
    ///   | 70 dB, 3/8    | 5     | 42 (40)           | -90.87 dB | 0.00895 dB | 92.02 dB     |
    ///   | 70 dB, 3/8    | 7     | 50 (50)           | -72.33 dB | 0.00966 dB | 70.41 dB     |
    ///   | 70 dB, 3/8    | 9     | 72 (64)           | -74.60 dB | 0.00929 dB | 73.46 dB     |
    ///   | 70 dB, 3/8    | 11    | 96 (86)           | -75.01 dB | 0.00975 dB | 72.48 dB     |
    ///   | 70 dB, 19/48  | 3     | 32 (32)           | -77.01 dB | 0.00816 dB | 73.97 dB     |
    ///   | 70 dB, 19/48  | 5     | 48 (48)           | -71.71 dB | 0.00983 dB | 70.10 dB     |
    ///   | 70 dB, 19/48  | 7     | 60 (58)           | -70.66 dB | 0.00964 dB | 70.36 dB     |
    ///   | 70 dB, 19/48  | 9     | 86 (76)           | -75.93 dB | 0.00973 dB | 81.04 dB     |
    ///   | 70 dB, 19/48  | 11    | 114 (102)         | -70.88 dB | 0.00998 dB | 70.01 dB     |
    ///
    /// The taps are for Q15 only: double, float and Q31 run k_guard_band_fir,
    /// which Q1.30 keeps (tests/test_iir_design.cpp). No row at 120 dB: Q1.14
    /// taps cannot carry it. Rounding each tap to Q1.14 (an LSB of 2^-14,
    /// -84.3 dB) adds an error response whose RMS over frequency is
    /// 2^-14 sqrt(N / 12), -76.5 dB at 72 taps, tens of dB above -120; the
    /// committed 120 dB rows, quantized so, attain -97.4 / -70.9 / -82.1 /
    /// -69.2 dB at orders 5 / 7 / 9 / 11 (pinned there). SampleRateTap's
    /// rational states Q15 at transparent the same way for its FIR stages (its
    /// PLAN.md section 6): format-limited, the 120 dB promise float's and
    /// Q31's.
    struct guard_band_fir_q15_design {
        int                     stop_db   = 0;
        int                     pass_num  = 0;
        int                     pass_den  = 1;
        int                     order     = 0;     ///< the Chebyshev II whose stop edge it matches
        double                  stop_edge = 0.0;   ///< fraction of its rate: the double row's
        std::span<const double> taps;              ///< the design, double; quantize it to Q1.14 as above
        double                  design_db   = 0.0; ///< the stopband the remez design targeted
        double                  ripple_db   = 0.0; ///< passband, peak to peak, attained by the Q1.14 taps
        double                  stopband_db = 0.0; ///< attained by the Q1.14 taps, from stop_edge

        constexpr int key() const noexcept { return order; }
    };

} // namespace tap::dsp

#include "tap/dsp/detail/iir_tables_data.h"
