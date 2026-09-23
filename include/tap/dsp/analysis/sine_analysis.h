/// @file sine_analysis.h
/// @brief Least-squares sine fitting for THD+N-style quality measurements.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Carried from SampleRateTap's test harness (tests/support/), promoted here
// as a shared measurement instrument: fit a sine of known (or tracked)
// frequency by least squares, subtract it, and score the residual. Used by
// converter quality suites (SampleRateTap, RatioTap) and available to any
// primitive that needs a single-tone SNR number.
//
// Sample types (Stage 3c of docs/audit-fft-and-code-smells.md): every entry
// point takes a contiguous range of Sample — std::span<const Sample>,
// std::vector<Sample>, a std::array — for Sample float, double, std::int16_t
// or std::int32_t. Floating samples are read as they are, exactly as before
// the instruments were typed (pinned bit for bit by test_analysis.cpp,
// `FloatingSpansAreBitIdenticalToThePreTemplateInstrument`); integer samples
// are read as Q0.15 / Q0.31 fractions of full scale through
// sample_traits<Sample>::k_sample_frac_bits, so a Q15 converter tail or a
// fixed-point FFT's time-domain output is scored in the same units (1.0 =
// full scale) as the floating profiles. The optional `exponent` is for data
// that carries a scale exponent, i.e. the output of a fixed-point real FFT
// (fft.h): the samples are read as x * 2^exponent, so a block-floating
// inverse's reconstruction is scored at its true level. The fit's numbers
// (amplitude, dc, residual_rms) are in those fraction units; snr_db is
// scale-free. The read is one exact power-of-two multiply per sample, so
// the hot loops stay the plain double loops they were.
#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <ranges>
#include <span>

#include "tap/dsp/sample_traits.h"

namespace tap::dsp::analysis {

    /// The sample types the instruments read: float and double as they are,
    /// std::int16_t and std::int32_t as Q0.15 / Q0.31 fractions of full scale.
    template <typename Sample>
    concept analysis_sample =
        std::floating_point<Sample> || std::same_as<Sample, std::int16_t> || std::same_as<Sample, std::int32_t>;

    /// A contiguous, sized range of analysis samples (std::span<const Sample>,
    /// std::vector<Sample>, std::array<Sample, N>, ...): what every instrument
    /// entry point takes, so a caller passes its buffer without spelling the
    /// span, and a std::span<const Sample> subrange still deduces.
    template <typename R>
    concept analysis_range = std::ranges::contiguous_range<R> && std::ranges::sized_range<R>
                             && analysis_sample<std::remove_cv_t<std::ranges::range_value_t<R>>>;

    namespace detail {

        /// Fraction bits of a sample type as the instruments read it: 0 for
        /// the floating types (the value is the fraction), the substrate's
        /// k_sample_frac_bits (15 / 31) for the fixed-point types.
        template <analysis_sample Sample>
        constexpr int fraction_bits() noexcept {
            if constexpr (std::floating_point<Sample>) {
                return 0;
            }
            else {
                return sample_traits<Sample>::k_sample_frac_bits;
            }
        }

        /// The multiplier that reads a sample as a fraction of full scale at
        /// scale exponent `exponent`: 2^(exponent - fraction bits). A power of
        /// two, so the product is exact (an int16 or int32 value fits a double
        /// mantissa), and for the floating types at exponent 0 it is 1.0, the
        /// identity: the floating path computes exactly what it did before.
        template <analysis_sample Sample>
        inline double sample_step(int exponent) noexcept {
            return std::ldexp(1.0, exponent - fraction_bits<Sample>());
        }

        /// The range as a std::span<const Sample>.
        template <analysis_range R>
        auto as_span(const R& x) noexcept {
            using sample = std::remove_cv_t<std::ranges::range_value_t<R>>;
            return std::span<const sample>(std::ranges::data(x), std::ranges::size(x));
        }

    } // namespace detail

    struct sine_fit {
        double amplitude    = 0.0;
        double phase        = 0.0;
        double dc           = 0.0;
        double residual_rms = 0.0;
        double freq_norm    = 0.0;
    };

    /// Fits x[i] ~ a*sin(w i) + b*cos(w i) + c by least squares (3x3 normal
    /// equations) at the known normalized frequency, then measures the residual.
    /// @param x         any contiguous range of float / double / Q15 / Q31 samples
    /// @param freq_norm the tone's frequency in cycles per sample
    /// @param exponent  scale exponent of the data (fixed-point FFT output): read as x * 2^exponent
    template <analysis_range R>
    sine_fit fit_sine(const R& x, double freq_norm, int exponent = 0) {
        const auto   xs      = detail::as_span(x);
        const double step    = detail::sample_step<typename decltype(xs)::value_type>(exponent);
        const double w       = 2.0 * std::numbers::pi * freq_norm;
        double       m[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        double       rhs[3]  = {0, 0, 0};
        for (std::size_t i = 0; i < xs.size(); ++i) {
            const double s        = std::sin(w * static_cast<double>(i));
            const double c        = std::cos(w * static_cast<double>(i));
            const double basis[3] = {s, c, 1.0};
            const double v        = static_cast<double>(xs[i]) * step;
            for (int r = 0; r < 3; ++r) {
                for (int q = 0; q < 3; ++q) {
                    m[r][q] += basis[r] * basis[q];
                }
                rhs[r] += basis[r] * v;
            }
        }
        // Gaussian elimination with partial pivoting.
        int order[3] = {0, 1, 2};
        for (int col = 0; col < 3; ++col) {
            int piv = col;
            for (int r = col + 1; r < 3; ++r) {
                if (std::abs(m[order[r]][col]) > std::abs(m[order[piv]][col])) {
                    piv = r;
                }
            }
            std::swap(order[col], order[piv]);
            const int p = order[col];
            for (int r = col + 1; r < 3; ++r) {
                const int    rr = order[r];
                const double f  = m[rr][col] / m[p][col];
                for (int q = col; q < 3; ++q) {
                    m[rr][q] -= f * m[p][q];
                }
                rhs[rr] -= f * rhs[p];
            }
        }
        double sol[3];
        for (int col = 2; col >= 0; --col) {
            const int p = order[col];
            double    v = rhs[p];
            for (int q = col + 1; q < 3; ++q) {
                v -= m[p][q] * sol[q];
            }
            sol[col] = v / m[p][col];
        }
        sine_fit fit;
        fit.amplitude = std::hypot(sol[0], sol[1]);
        fit.phase     = std::atan2(sol[1], sol[0]);
        fit.dc        = sol[2];
        double sq     = 0.0;
        for (std::size_t i = 0; i < xs.size(); ++i) {
            const double s = std::sin(w * static_cast<double>(i));
            const double c = std::cos(w * static_cast<double>(i));
            const double r = static_cast<double>(xs[i]) * step - (sol[0] * s + sol[1] * c + sol[2]);
            sq += r * r;
        }
        fit.residual_rms = std::sqrt(sq / static_cast<double>(xs.size()));
        fit.freq_norm    = freq_norm;
        return fit;
    }

    /// Like fitSine, but refines the frequency first (a few iterations comparing
    /// the fitted phase of the two window halves). A converter's rate estimate
    /// converges asymptotically, so the tail of a run can sit a fraction of a ppm
    /// off the nominal ratio; a rigid fixed-frequency fit would book that
    /// (inaudible) offset as residual. Tracking the fundamental is standard
    /// THD-analyzer practice.
    /// @param x, exponent  as for fit_sine
    /// @param freq_norm_guess  the nominal frequency in cycles per sample; the fit refines it
    template <analysis_range R>
    sine_fit fit_sine_tracked(const R& x, double freq_norm_guess, int exponent = 0) {
        const auto        xs   = detail::as_span(x);
        double            f    = freq_norm_guess;
        const std::size_t half = xs.size() / 2;
        for (int iter = 0; iter < 4; ++iter) {
            const sine_fit a = fit_sine(xs.first(half), f, exponent);
            const sine_fit b = fit_sine(xs.subspan(half), f, exponent);
            // b.phase is relative to the second half's start; predict it from a.
            const double two_pi    = 2.0 * std::numbers::pi;
            const double predicted = a.phase + two_pi * f * static_cast<double>(half);
            const double dphi      = std::remainder(b.phase - predicted, two_pi);
            f += dphi / (two_pi * static_cast<double>(half));
        }
        return fit_sine(xs, f, exponent);
    }

    /// Signal-to-(residual) ratio in dB for a fitted sine.
    inline double snr_db(const sine_fit& f) {
        return 10.0 * std::log10((f.amplitude * f.amplitude * 0.5) / (f.residual_rms * f.residual_rms));
    }

} // namespace tap::dsp::analysis
