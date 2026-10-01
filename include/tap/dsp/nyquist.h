/// @file nyquist.h
/// @brief L-th-band (Nyquist) FIR design: the lowpass whose every L-th tap is zero.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The design math of the SampleRateTap family's `rational` engine (its
// PLAN.md, section 2.2), landed here first per the substrate rule: a stage
// that changes the rate by a small factor L (or 1/L) is one L-th-band filter,
// and decimate.h may adopt the same designer later without touching the
// engine. Published literature only: the L-th-band structure is Mintzer 1982
// ("On half-band, third-band and Nth-band FIR filters") and Vaidyanathan,
// Multirate Systems and Filter Banks, section 4.6; the window and its fit
// are kaiser.h's (Kaiser 1974, harris).
//
// Contract, as numbers:
//   - Length N = 2mL - 1 for an integer m >= 1 (nyquist_length), centre
//     c = mL - 1 (nyquist_centre), the design on the HIGHER rate of the
//     stage: the filter runs at F_hi = L * r where r is the lower rate.
//   - h[c] == 1.0 exactly and h[c +- kL] == 0.0 exactly for k = 1 .. m - 1:
//     the Kaiser-windowed sinc with cutoff exactly pi/L has its zeros at
//     n - c = +-L, +-2L, ... by construction (the window multiplies them, so
//     they stay zero); the designer writes them as exact zeros rather than the
//     ~1e-16 libm leaves, because the contract says 0 and the fixed-point
//     tables are bit-pinned. Nonzero taps: 2m(L - 1) + 1
//     (nyquist_nonzero_taps). N is odd, so the group delay (N - 1) / 2 is an
//     integer number of samples at F_hi.
//   - Normalization: every polyphase branch (the taps i == j mod L) has DC
//     gain 1 within double rounding, and the branch holding the centre is
//     the single tap h[c] == 1.0, exactly: an interpolator's branch 0 is a
//     copy. The whole filter therefore sums to L, kaiser.h's
//     design_prototype convention. A decimator by L uses h / L (sum 1); the
//     fixed-point rows land on exact unity through quantize.h.
//   - Symmetric transition: the response is antisymmetric about F_hi / (2L)
//     = r / 2, so a passband edge f_p and a stopband edge f_s satisfy
//     f_p + f_s = r (the Nyquist property in frequency, sum over k of
//     H(w - 2 pi k / L) == 1, nyquist_response_db pins it). The engine's
//     coverage rule asks for f_s <= r - f_p; the structure gives equality.
//   - The spec search (search_nyquist_m) steps m upward until the worst
//     stopband response on a fine grid is <= -(A + margin) dB, the bridge
//     engine's criterion, so tap counts step by 2L. Direct DFT by rotator:
//     no FFT, no dependency; design time only.
//
// Runtime double, like every designer here (kaiser.h's design note); the
// designer itself is noexcept and allocation-free, the search allocates.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>

#include "tap/dsp/kaiser.h"

namespace tap::dsp {

    /// N = 2mL - 1: the length of the L-th-band filter with m taps per branch.
    constexpr std::size_t nyquist_length(std::size_t num_phases, std::size_t m) noexcept {
        return 2 * m * num_phases - 1;
    }

    /// The centre tap's index, mL - 1 == (N - 1) / 2.
    constexpr std::size_t nyquist_centre(std::size_t num_phases, std::size_t m) noexcept {
        return m * num_phases - 1;
    }

    /// The taps that are not structurally zero: 2m(L - 1) + 1.
    constexpr std::size_t nyquist_nonzero_taps(std::size_t num_phases, std::size_t m) noexcept {
        return 2 * m * (num_phases - 1) + 1;
    }

    /// True when n is a valid L-th-band length, 2mL - 1 for some m >= 1.
    constexpr bool is_nyquist_length(std::size_t num_phases, std::size_t n) noexcept {
        return num_phases >= 1 && n >= 2 * num_phases - 1 && (n + 1) % (2 * num_phases) == 0;
    }

    /// Designs the L-th-band Kaiser-windowed sinc lowpass into h.
    ///
    /// @param h          output, size nyquist_length(num_phases, m) for some m >= 1
    ///                   (is_nyquist_length); any other size is a contract
    ///                   violation and h is left untouched
    /// @param num_phases L: the band index and the polyphase count
    /// @param beta       Kaiser shape parameter (kaiser_beta)
    ///
    /// Centre tap exactly 1.0, every L-th tap from it exactly 0.0, each
    /// branch normalized to DC gain 1 (sum(h) == L); see the file header.
    inline void design_nyquist(std::span<double> h, std::size_t num_phases, double beta) noexcept {
        const std::size_t n = h.size();
        if (!is_nyquist_length(num_phases, n)) {
            return;
        }
        const std::size_t c = (n - 1) / 2;
        if (c == 0) { // L == 1, m == 1: the one-tap identity
            h[0] = 1.0;
            return;
        }
        const double i0_beta = bessel_i0(beta);
        // Kaiser window about the centre, half the Bessel series through the
        // mirror (the argument of tap n-1-i is the exact negation of tap i's).
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t mirror = n - 1 - i;
            if (mirror < i) {
                h[i] = h[mirror];
                continue;
            }
            const double u = (static_cast<double>(i) - static_cast<double>(c)) / static_cast<double>(c);
            h[i]           = bessel_i0(beta * std::sqrt(std::max(0.0, 1.0 - u * u)));
        }
        // Windowed sinc with cutoff pi/L: (1/L) sinc((i - c)/L). The structural
        // zeros (i == c mod L, i != c) and the centre are written exactly.
        const double inv_l = 1.0 / static_cast<double>(num_phases);
        for (std::size_t i = 0; i < n; ++i) {
            if (i == c) {
                h[i] = 1.0;
            }
            else if ((i + num_phases - (c % num_phases)) % num_phases == 0) {
                h[i] = 0.0;
            }
            else {
                const double t = (static_cast<double>(i) - static_cast<double>(c)) * inv_l;
                h[i]           = inv_l * sinc(t) * (h[i] / i0_beta);
            }
        }
        // Per-branch DC normalization. The centre's branch is the single tap
        // h[c] == 1.0 already; every other branch is scaled to sum 1.
        for (std::size_t j = 0; j < num_phases; ++j) {
            if (j == c % num_phases) {
                continue;
            }
            double sum = 0.0;
            for (std::size_t i = j; i < n; i += num_phases) {
                sum += h[i];
            }
            const double gain = 1.0 / sum;
            for (std::size_t i = j; i < n; i += num_phases) {
                h[i] *= gain;
            }
        }
    }

    /// Magnitude response of an L-th-band design in dB, 0 dB at DC, at the
    /// normalized frequency f_norm = f / F_hi (F_hi = L * r, the filter's own
    /// rate). Direct DFT by rotator; design time only.
    inline double nyquist_response_db(std::span<const double> h, std::size_t num_phases, double f_norm) noexcept {
        const double th   = -2.0 * std::numbers::pi * f_norm;
        const double th_c = std::cos(th);
        const double th_s = std::sin(th);
        double       rc = 1.0, rs = 0.0; // e^{j th i}, advanced per tap
        double       re = 0.0, im = 0.0;
        for (const double v : h) {
            re += v * rc;
            im += v * rs;
            const double nrc = rc * th_c - rs * th_s;
            rs               = rs * th_c + rc * th_s;
            rc               = nrc;
        }
        const double mag = std::sqrt(re * re + im * im) / static_cast<double>(num_phases);
        return 20.0 * std::log10(std::max(mag, 1e-300));
    }

    /// Worst stopband response in dB of an L-th-band design, swept on a grid
    /// from the stopband edge to F_hi / 2. The passband edge is given as a
    /// fraction p of the LOWER rate r (f_p = p * r, 0 < p < 1/2); the stopband
    /// edge is then r - f_p, i.e. (1 - p) / L of F_hi. Design time only.
    inline double nyquist_worst_stopband_db(std::span<const double> h, std::size_t num_phases, double passband_frac,
                                            std::size_t grid_points = 1024) noexcept {
        const double f_stop = (1.0 - passband_frac) / static_cast<double>(num_phases);
        double       worst  = -1e9;
        for (std::size_t g = 0; g < grid_points; ++g) {
            const double f = f_stop + (0.5 - f_stop) * static_cast<double>(g) / static_cast<double>(grid_points - 1);
            worst          = std::max(worst, nyquist_response_db(h, num_phases, f));
        }
        return worst;
    }

    /// The smallest m >= 1 whose L-th-band design meets the stopband spec with
    /// the margin on nyquist_worst_stopband_db's grid: the bridge engine's
    /// ">= 1 dB margin on a fine grid" criterion. Returns 0 when no m up to
    /// max_m does. Allocates (one candidate at a time); design time only.
    inline std::size_t search_nyquist_m(std::size_t num_phases, double passband_frac, double stopband_atten_db,
                                        double margin_db = 1.0, std::size_t max_m = 256) {
        const double        beta = kaiser_beta(stopband_atten_db);
        std::vector<double> h;
        for (std::size_t m = 1; m <= max_m; ++m) {
            h.assign(nyquist_length(num_phases, m), 0.0);
            design_nyquist(h, num_phases, beta);
            if (nyquist_worst_stopband_db(h, num_phases, passband_frac) <= -(stopband_atten_db + margin_db)) {
                return m;
            }
        }
        return 0;
    }

} // namespace tap::dsp
