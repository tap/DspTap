/// @file biquad.h
/// @brief The second-order-section (biquad) cascade kernel: transposed direct form II over float and double,
/// its group delay, and its exact-zero flush window.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The recursive substrate SampleRateTap's IIR plan places in DspTap (its I1:
// "the biquad / first-order cascade, per sample_traits format"), landed first
// under butterworth.h. A cascade of S sections, each
//
//     H_s(z) = (b0 + b1 z^-1 + b2 z^-2) / (1 + a1 z^-1 + a2 z^-2),
//
// run in transposed direct form II (Jackson, "Digital Filters and Signal
// Processing", 3rd ed., 1996, sec. 11.4; the form the plan's I8 names for
// float):
//
//     y = b0 x + s1;   s1 <- b1 x - a1 y + s2;   s2 <- b2 x - a2 y.
//
// Contract, as numbers:
//   - Formats: double (the golden model) and float. Float keeps its
//     coefficients, products and states in float (the plan's I14: a biquad
//     is not the fragile kind of recursion that DspTap otherwise runs in
//     double; IIR-0 measured float states 28 times cheaper than double on
//     the soft-double M33 and M4F). Q15 / Q31 are not offered yet: the
//     fixed-point design (direct form I, Q2.30 coefficients, 64-bit
//     accumulation, magnitude truncation, L-infinity gain distribution) is
//     the plan's I8 and lands with it; a static_assert names it.
//   - Frames are interleaved, Ch channels (a compile-time count); each
//     channel runs its own state through the same coefficients.
//   - Section order is the caller's; a designer orders by increasing pole
//     radius (butterworth.h), the order the plan's I8 requires for fixed
//     point and harmless in float.
//   - Group delay (group_delay_samples): the sum over sections of the
//     numerator's minus the denominator's, Re(sum k p_k e^{-iwk} /
//     sum p_k e^{-iwk}), evaluated in double from the double design, for
//     0 <= f < 1/2 cycles per sample. A numerator with a zero on the unit
//     circle at f is singular there (a lowpass's Nyquist zero: f < 1/2).
//   - Non-finite input poisons the state until reset() (the plan's I9: an
//     honest limit, no per-sample check).
//
// The tail: two policies, a template argument (iir_tail).
//
// Why there are two (measured, x86-64 and by construction): under the
// host's flush-to-zero a float biquad does NOT in general decay to exact
// zero. FTZ flushes each product on its own, so near the smallest normal
// value a section can lose a2 y while s1 survives; with |a1| > 1 (any pole
// pair with a real part above 1/2) s1 then grows by |a1| until a2 y is
// normal again, and the state settles into a limit cycle a few times the
// smallest normal value. An order-4 Butterworth at fc = 0.1 holds one
// forever in transposed DF II, DF I and (at other cutoffs) DF II, float and
// double alike; the plan's IIR-0 measured exact decay for its first-order
// allpass sections only, whose products can only shrink when flushed. So
// the plan's I6 contract (after window_frames() zeros the state is exactly
// zero) needs the kernel's help for biquads:
//
//   - iir_tail::exact_zero (the default): after each section's update, if
//     both its state words are below T in magnitude, both are set to zero
//     (branch-free: a compare and two selects). T = 2^(digits + 2) times
//     the smallest normal value: 2^-100 in float (-602 dB re full scale),
//     2^-967 in double. Every product of a state above T with a
//     coefficient stays normal, so FTZ never acts on a live state, and the
//     joint test cannot strand one word (zeroing a single word is what
//     FTZ does, and what cycles). The tail reaches exact zero within
//     flush_window_frames<S>(sections) zero frames with or without the
//     host's FTZ, so this policy states no subnormal precondition. Below
//     T the filter is not linear: a state that small is dropped. Cost:
//     about 26 % more time per sample on x86-64 at order 4, mono, float
//     (one dependent chain; the instruction counts on CI's legs are in
//     bench/README.md).
//   - iir_tail::free_running: the plain recursion. Under FTZ the tail can
//     hold a limit cycle at a few times the smallest normal value, and
//     without it a few to some hundred multiples of the smallest subnormal
//     (measured up to 114 at order 8, fc = 0.005, float), so the stage
//     has no exact window: a stage built on it offers no window_frames()
//     and a chain holding it cannot flush (chain.h's static_assert). Its
//     process() states the plan's I9 precondition, subnormals_flushed<S>(),
//     because without FTZ x86 spends its silence in subnormal arithmetic.
//
// The flush window (exact_zero; the plan's I6: the count of zero input
// frames after which the state is exactly zero, so chain.h's flush
// contract holds unchanged). With zero input the cascade is x' = A x, A
// the 2S x 2S state-transition matrix. The state holds at most
// max_finite / G (the state headroom), G = sup_m ||A^m||_inf. The window
// is the smallest n with
//
//     ln(max_finite) + ln ||A^n||_inf  <  ln(T),
//
// because then every later state obeys ||x_m|| <= ||A^n|| ||A^(m-n)|| ||x_0||
// <= ||A^n|| G max_finite / G < T, where every section is zeroed; plus a
// margin of one frame per section (the zeroed sections pass silence down
// the cascade one section per frame at most; IIR-0 found the same margin).
// n is found by binary lifting over A^(2^j), kept normalized with a
// separate log scale, so the design costs O(log n) matrix products even at
// a cutoff of a few hertz in double (n in the millions). The count is a
// design-time constant.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include "tap/dsp/subnormals.h"

namespace tap::dsp {

    /// What a recursive kernel does with its tail (the file header).
    enum class iir_tail {
        exact_zero,   ///< zero a section's state jointly below T: an exact flush window
        free_running, ///< the plain recursion: no exact window, the host's FTZ required
    };

    /// T, the exact_zero policy's threshold: 2^(digits + 2) times the smallest normal value.
    template <typename S>
    constexpr S biquad_tail_threshold() noexcept {
        S t = std::numeric_limits<S>::min();
        for (int i = 0; i < std::numeric_limits<S>::digits + 2; ++i) {
            t *= S(2);
        }
        return t;
    }

    /// One second-order section, normalized so a0 = 1.
    struct biquad_coefficients {
        double b0 = 1.0, b1 = 0.0, b2 = 0.0;
        double a1 = 0.0, a2 = 0.0;
    };

    /// Group delay of a section cascade in samples at f cycles per sample, 0 <= f < 1/2.
    template <std::size_t S>
    double group_delay_samples(const std::array<biquad_coefficients, S>& sections, double f) noexcept {
        TAP_EXPECTS(f >= 0.0 && f < 0.5);
        constexpr double           pi = 3.14159265358979323846;
        const std::complex<double> z1 = std::polar(1.0, -2.0 * pi * f);
        const std::complex<double> z2 = z1 * z1;
        double                     gd = 0.0;
        for (const auto& c : sections) {
            const std::complex<double> num  = c.b0 + c.b1 * z1 + c.b2 * z2;
            const std::complex<double> dnum = c.b1 * z1 + 2.0 * c.b2 * z2;
            const std::complex<double> den  = 1.0 + c.a1 * z1 + c.a2 * z2;
            const std::complex<double> dden = c.a1 * z1 + 2.0 * c.a2 * z2;
            gd += (dnum / num).real() - (dden / den).real();
        }
        return gd;
    }

    namespace detail {

        /// A square matrix kept as M * e^log_scale, M normalized to ||M||_inf = 1.
        template <std::size_t D>
        struct scaled_matrix {
            std::array<double, D * D> m{};
            double                    log_scale = 0.0;

            static double norm_inf(const std::array<double, D * D>& a) noexcept {
                double n = 0.0;
                for (std::size_t r = 0; r < D; ++r) {
                    double row = 0.0;
                    for (std::size_t c = 0; c < D; ++c) {
                        row += std::fabs(a[r * D + c]);
                    }
                    n = std::max(n, row);
                }
                return n;
            }
            void normalize() noexcept {
                const double n = norm_inf(m);
                if (n > 0.0) {
                    for (double& v : m) {
                        v /= n;
                    }
                    log_scale += std::log(n);
                }
                else {
                    log_scale = -std::numeric_limits<double>::infinity();
                }
            }
            /// ln ||this||_inf.
            double log_norm() const noexcept { return log_scale; }

            friend scaled_matrix operator*(const scaled_matrix& a, const scaled_matrix& b) noexcept {
                scaled_matrix r;
                for (std::size_t i = 0; i < D; ++i) {
                    for (std::size_t k = 0; k < D; ++k) {
                        const double aik = a.m[i * D + k];
                        for (std::size_t j = 0; j < D; ++j) {
                            r.m[i * D + j] += aik * b.m[k * D + j];
                        }
                    }
                }
                r.log_scale = a.log_scale + b.log_scale;
                r.normalize();
                return r;
            }
        };

        /// The zero-input state-transition matrix of a TDF-II cascade; state
        /// order (s1, s2) per section, sections in processing order.
        template <std::size_t S>
        std::array<double, 4 * S * S> cascade_transition(const std::array<biquad_coefficients, S>& sec) noexcept {
            constexpr std::size_t     D = 2 * S;
            std::array<double, D * D> a{};
            // x_s: the input to section s, as a row vector over the old state.
            std::array<double, D> x{}; // section 0's input is the zero input
            for (std::size_t s = 0; s < S; ++s) {
                const auto& c = sec[s];
                // y = b0 x + s1 (old)
                std::array<double, D> y = x;
                for (double& v : y) {
                    v *= c.b0;
                }
                y[2 * s] += 1.0;
                // s1' = b1 x - a1 y + s2 (old); s2' = b2 x - a2 y
                for (std::size_t j = 0; j < D; ++j) {
                    a[(2 * s) * D + j]     = c.b1 * x[j] - c.a1 * y[j];
                    a[(2 * s + 1) * D + j] = c.b2 * x[j] - c.a2 * y[j];
                }
                a[(2 * s) * D + 2 * s + 1] += 1.0;
                x = y; // the next section's input is this section's output
            }
            return a;
        }

    } // namespace detail

    /// Zero input frames after which an iir_tail::exact_zero cascade's state
    /// is exactly zero, from any state format S can hold (the file header).
    /// Design-time: O(log n) products of 2S x 2S matrices.
    template <typename S, std::size_t N>
    std::size_t flush_window_frames(const std::array<biquad_coefficients, N>& sections) {
        constexpr std::size_t D = 2 * N;
        using mat               = detail::scaled_matrix<D>;
        mat a;
        a.m = detail::cascade_transition(sections);
        a.normalize();

        const double start  = std::log(static_cast<double>(std::numeric_limits<S>::max()));
        const double target = std::log(static_cast<double>(biquad_tail_threshold<S>()));
        auto         done   = [&](const mat& p) { return start + p.log_norm() < target; };

        // Powers A^(2^j) until one alone satisfies the condition.
        std::vector<mat> pow2{a};
        while (!done(pow2.back()) && pow2.size() < 63) {
            pow2.push_back(pow2.back() * pow2.back());
        }
        // Binary lifting: the largest n built from those powers that is not yet done.
        mat         acc;
        bool        have = false;
        std::size_t n    = 0;
        for (std::size_t j = pow2.size(); j-- > 0;) {
            const mat cand = have ? acc * pow2[j] : pow2[j];
            if (!done(cand)) {
                acc  = cand;
                have = true;
                n += std::size_t{1} << j;
            }
        }
        // Step forward until n satisfies the condition (the lifting found
        // the last n built from powers of two that does not).
        mat p = have ? acc * a : a;
        ++n;
        while (!done(p)) {
            p = p * a;
            ++n;
        }
        return n + N; // the margin: one frame per section
    }

    /// The cascade kernel: N sections, Ch interleaved channels, float or double.
    template <typename S, std::size_t N, std::size_t Ch, iir_tail Tail = iir_tail::exact_zero>
    class basic_biquad_cascade {
        static_assert(std::is_same_v<S, float> || std::is_same_v<S, double>,
                      "basic_biquad_cascade: float and double; Q15 / Q31 are the IIR plan's I8 (fixed-point "
                      "direct form I), not implemented yet");
        static_assert(N >= 1 && Ch >= 1);

      public:
        using sample = S;

        basic_biquad_cascade() = default;
        explicit basic_biquad_cascade(const std::array<biquad_coefficients, N>& sec) noexcept { set(sec); }

        void set(const std::array<biquad_coefficients, N>& sec) noexcept {
            for (std::size_t s = 0; s < N; ++s) {
                m_c[s] = {static_cast<S>(sec[s].b0), static_cast<S>(sec[s].b1), static_cast<S>(sec[s].b2),
                          static_cast<S>(sec[s].a1), static_cast<S>(sec[s].a2)};
            }
        }

        /// n interleaved frames; in and out may be the same buffer.
        void process(const S* in, std::size_t n, S* out) noexcept {
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t ch = 0; ch < Ch; ++ch) {
                    S x = in[i * Ch + ch];
                    for (std::size_t s = 0; s < N; ++s) {
                        const coeffs& c  = m_c[s];
                        S*            st = m_state[s][ch].data();
                        const S       y  = c.b0 * x + st[0];
                        S             s0 = c.b1 * x - c.a1 * y + st[1];
                        S             s1 = c.b2 * x - c.a2 * y;
                        if constexpr (Tail == iir_tail::exact_zero) {
                            const bool below = (s0 < k_t && s0 > -k_t) && (s1 < k_t && s1 > -k_t);
                            s0               = below ? S(0) : s0;
                            s1               = below ? S(0) : s1;
                        }
                        st[0] = s0;
                        st[1] = s1;
                        x     = y;
                    }
                    out[i * Ch + ch] = x;
                }
            }
        }

        void reset() noexcept {
            for (auto& sec : m_state) {
                for (auto& st : sec) {
                    st = {S(0), S(0)};
                }
            }
        }

        /// True when every state word is exactly zero (silence from here on).
        bool is_silent() const noexcept {
            for (const auto& sec : m_state) {
                for (const auto& st : sec) {
                    if (st[0] != S(0) || st[1] != S(0)) {
                        return false;
                    }
                }
            }
            return true;
        }

      private:
        static constexpr S k_t = biquad_tail_threshold<S>();
        struct coeffs {
            S b0, b1, b2, a1, a2;
        };
        std::array<coeffs, N>                           m_c{};
        std::array<std::array<std::array<S, 2>, Ch>, N> m_state{};
    };

} // namespace tap::dsp
