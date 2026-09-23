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
// or std::int32_t (exactly those four; long double is not admitted, since it
// would be read through a narrowing cast). Floating samples are read as they
// are, exactly as before the instruments were typed: the float and double
// instantiations are bit-identical to the pre-template std::span<const float>
// instrument under -ffp-contract=off (the flag test_analysis_typed.cpp's
// target carries, `FloatingSpansAreBitIdenticalToThePreTemplateInstrument`)
// and on compilers that contract per statement (clang at `on` and `fast`,
// MSVC); under GCC's default cross-statement -ffp-contract=fast on FMA
// hardware they differ from it by at most one ulp (measured 2026-09-23,
// g++ 13.3 -O3 -mfma: amplitude and phase of the tracked float fit each moved
// one ulp; g++ with -ffp-contract=off, g++ without FMA, and clang++ -mfma at
// `on` and `fast` were identical). Integer samples are read as Q0.15 / Q0.31
// fractions of full scale through sample_traits<Sample>::k_sample_frac_bits,
// so a Q15 converter tail or a fixed-point FFT's time-domain output is scored
// in the same units (1.0 = full scale) as the floating profiles. The optional
// `exponent` is for data that carries a scale exponent, i.e. the output of a
// fixed-point real FFT (fft.h): the fit is reported as if the samples were
// x * 2^exponent, so a block-floating inverse's reconstruction is scored at
// its true level. The fit's numbers (amplitude, dc, residual_rms) are in
// those fraction units; snr_db is scale-free.
//
// Source compatibility: every call site that passed a std::span<const float>,
// a vector or an array by name still compiles unchanged. The one shape that
// no longer deduces is a braced span built in the call, `fit_sine({ptr, n},
// nu)`; spell it `fit_sine(std::span<const float>(ptr, n), nu)`. No caller in
// DspTap or MuTap used it.
//
// How the read is done matters for bit identity, so it is stated: a
// floating sample enters the arithmetic as the bare static_cast<double> it
// always was (no multiply, not even by 1.0: under fp-contraction a product
// feeding a subtraction fuses differently from a plain operand, and the
// pre-template bits would move by an ulp on hosts that contract, which the
// macOS arm64 CI leg showed); an integer sample enters as its exact fraction
// (an int16 / int32 times 2^-15 / 2^-31 is exact, so fusing that product
// into a later add cannot change a bit either); and the exponent is applied
// to the three magnitude fields after the fit, an exact power-of-two scale,
// never to the samples. The hot loops stay the plain double loops they were.
//
// Includes: the range concept is written over std::data / std::size, which
// <span> already provides ([iterator.range]), rather than over <ranges> or
// <iterator>: measured on one TU including only this header (-std=c++20 -O2
// -c, best of three, 2026-09-23), g++ 13.3 takes 0.19 s on the pre-Stage-3c
// header, 0.40 s with <ranges>, 0.35 s with <iterator> and 0.19 s as written;
// clang++ 18.1 0.23 / 0.58 / 0.44 / 0.23 s.
#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <span>
#include <type_traits>
#include <utility>

#include "tap/dsp/sample_traits.h"

namespace tap::dsp::analysis {

    /// The sample types the instruments read: float and double as they are,
    /// std::int16_t and std::int32_t as Q0.15 / Q0.31 fractions of full scale.
    /// Exactly these four: long double is not admitted (it would be narrowed).
    template <typename Sample>
    concept analysis_sample = std::same_as<Sample, float> || std::same_as<Sample, double>
                              || std::same_as<Sample, std::int16_t> || std::same_as<Sample, std::int32_t>;

    /// The element type a contiguous range holds, from what std::data(r) points at.
    template <typename R>
    using range_sample_t = std::remove_cv_t<std::remove_pointer_t<decltype(std::data(std::declval<const R&>()))>>;

    /// A contiguous, sized range of analysis samples (std::span<const Sample>,
    /// std::vector<Sample>, std::array<Sample, N>, a C array, ...): what every
    /// instrument entry point takes, so a caller passes its buffer without
    /// spelling the span, and a std::span<const Sample> subrange still
    /// deduces. Spelled over std::data / std::size, which <span> provides
    /// ([iterator.range]): a type whose std::data yields a pointer to its
    /// elements and whose std::size yields their count is contiguous by that
    /// contract, so the header pulls neither <ranges> nor <iterator> (whose
    /// std::contiguous_iterator would restate it); each costs more to compile
    /// than the rest of this header (file comment).
    template <typename R>
    concept analysis_range = requires(const R& r) {
        { std::data(r) } -> std::convertible_to<const range_sample_t<R>*>;
        { std::size(r) } -> std::convertible_to<std::size_t>;
    } && analysis_sample<range_sample_t<R>>;

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

        /// One sample as a double in the instrument's units: a floating sample
        /// as it is, the bare cast the pre-template instrument evaluated; an
        /// integer sample as its exact Q0.15 / Q0.31 fraction (a 16- or 32-bit
        /// integer times a power of two is exact in a double).
        template <analysis_sample Sample>
        constexpr double to_fraction(Sample v) noexcept {
            if constexpr (std::floating_point<Sample>) {
                return static_cast<double>(v);
            }
            else {
                constexpr double unit = 1.0 / static_cast<double>(std::int64_t{1} << fraction_bits<Sample>());
                return static_cast<double>(v) * unit;
            }
        }

        /// Reports a fit made on x as the fit of x * 2^exponent: the three
        /// magnitude fields scale exactly (a power of two), the phase and the
        /// frequency do not change. A no-op at exponent 0, so the default path
        /// touches nothing.
        template <typename Fit>
        inline void apply_exponent(Fit& fit, int exponent) noexcept {
            if (exponent != 0) {
                const double k = std::ldexp(1.0, exponent);
                fit.amplitude *= k;
                fit.dc *= k;
                fit.residual_rms *= k;
            }
        }

        /// The range as a std::span<const Sample>.
        template <analysis_range R>
        auto as_span(const R& x) noexcept {
            return std::span<const range_sample_t<R>>(std::data(x), static_cast<std::size_t>(std::size(x)));
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
    /// @param exponent  scale exponent of the data (fixed-point FFT output): the fit is reported for x * 2^exponent
    template <analysis_range R>
    sine_fit fit_sine(const R& x, double freq_norm, int exponent = 0) {
        const auto   xs      = detail::as_span(x);
        const double w       = 2.0 * std::numbers::pi * freq_norm;
        double       m[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        double       rhs[3]  = {0, 0, 0};
        for (std::size_t i = 0; i < xs.size(); ++i) {
            const double s        = std::sin(w * static_cast<double>(i));
            const double c        = std::cos(w * static_cast<double>(i));
            const double basis[3] = {s, c, 1.0};
            for (int r = 0; r < 3; ++r) {
                for (int q = 0; q < 3; ++q) {
                    m[r][q] += basis[r] * basis[q];
                }
                rhs[r] += basis[r] * detail::to_fraction(xs[i]);
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
            const double r = detail::to_fraction(xs[i]) - (sol[0] * s + sol[1] * c + sol[2]);
            sq += r * r;
        }
        fit.residual_rms = std::sqrt(sq / static_cast<double>(xs.size()));
        fit.freq_norm    = freq_norm;
        detail::apply_exponent(fit, exponent);
        return fit;
    }

    /// Like fitSine, but refines the frequency first (a few iterations comparing
    /// the fitted phase of the two window halves). A converter's rate estimate
    /// converges asymptotically, so the tail of a run can sit a fraction of a ppm
    /// off the nominal ratio; a rigid fixed-frequency fit would book that
    /// (inaudible) offset as residual. Tracking the fundamental is standard
    /// THD-analyzer practice.
    /// @param x                as for fit_sine
    /// @param freq_norm_guess  the nominal frequency in cycles per sample; the fit refines it
    /// @param exponent         as for fit_sine
    template <analysis_range R>
    sine_fit fit_sine_tracked(const R& x, double freq_norm_guess, int exponent = 0) {
        const auto        xs   = detail::as_span(x);
        double            f    = freq_norm_guess;
        const std::size_t half = xs.size() / 2;
        for (int iter = 0; iter < 4; ++iter) {
            const sine_fit a = fit_sine(xs.first(half), f); // the phases do not depend on the exponent
            const sine_fit b = fit_sine(xs.subspan(half), f);
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
