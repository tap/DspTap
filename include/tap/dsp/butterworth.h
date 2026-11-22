/// @file butterworth.h
/// @brief Butterworth low-pass: the even-order design (bilinear, pre-warped) and a 1:1 chain stage over it.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// A recursive 1:1 chain stage over the IIR substrate (iir_kernels.h's biquad
// cascade, iir_design.h's sections and group delay, chain.h's nominal
// latency), under SampleRateTap's IIR plan's contracts for a recursive stage
// in a chain: latency I5 and flush I6. Subnormals are the host's (I9,
// subnormals.h): the host sets flush-to-zero; the kernel only reads it, as
// a debug-build precondition.
//
// Design (Oppenheim & Schafer, "Discrete-Time Signal Processing", 3rd ed.,
// 2010, sec. 7.2: the bilinear transform with frequency pre-warping, and
// the Butterworth prototype's pole placement). For order N (even) and cutoff
// fc in cycles per sample (the -3 dB point), K = tan(pi fc); the analogue
// prototype's poles pair into N/2 sections s^2 + alpha_k s + 1 with
// alpha_k = 2 sin((2k + 1) pi / (2N)), and the bilinear map gives each
//
//     b = K^2 [1, 2, 1] / a0,   a1 = 2 (K^2 - 1) / a0,   a2 = (1 - alpha_k K + K^2) / a0,
//     a0 = 1 + alpha_k K + K^2.
//
// Each section has unity gain at DC and a double zero at Nyquist. The
// sections run in increasing pole radius (lowest Q first), the order the
// plan's I8 requires for fixed point. The design is runtime double at
// construction, like every DspTap designer (kaiser.h, nyquist.h): a few
// coefficients through libm's tan and sin. The plan's constexpr designers
// (I4) are for its named tables; a Butterworth with a run-time cutoff has
// none.
//
// Contract, as numbers:
//   - Order N even, 2 ... 16; cutoff 0 < fc < 1/2 cycles per sample
//     (butterworth_supports_cutoff; the constructor rejects anything else
//     through raise<std::invalid_argument>). Magnitude |H(fc)| = 1/sqrt(2)
//     (-3.0103 dB), |H(0)| = 1, H(1/2) = 0; monotone in between.
//   - A 1:1 sync_stage (chain.h): k_up = k_down = 1, channels() = Ch, any
//     Ch >= 1, frames interleaved; process() writes exactly n frames
//     (outputs_for(n) = n), noexcept, allocation-free, bit-identical for any
//     chunking (the state is per channel and per sample); across channel
//     counts and chain.h's flush and process paths, given fft.h's
//     fp-contraction policy (contraction off; fft_fir.h says why). The
//     arithmetic is iir_kernels.h's floating biquad cascade, transposed
//     direct form II, one state pair per section and channel.
//   - Latency (the plan's I5, chain.h's reports_nominal_latency): a
//     recursive stage has no single delay, so it never provides
//     latency_output_frames(). It reports three numbers, in frames (x 1 /
//     rate for seconds):
//       group_delay_output_frames(f): -d(phase)/d(omega) at f cycles per
//         sample, 0 <= f < 1/2 (iir_design.h's sections_group_delay, double);
//       nominal_latency_output_frames(): the group delay at DC, the delay a
//         caller gives a parallel path (a dry signal, a sibling channel) to
//         line it up with this one;
//       latency_deviation_output_frames(): the passband group-delay
//         deviation, the largest |group_delay(f) - nominal| over 0 <= f <=
//         fc, on a 1,024-point grid. It is the alignment error a caller
//         accepts by compensating the nominal alone: a passband component
//         still arrives up to this many frames off the delayed parallel
//         path, so a sum of the two cancels where that offset is half the
//         component's period. A Butterworth's group delay rises toward fc,
//         so the deviation is its excess at the band edge (order 4: 2.51
//         frames on a nominal of 4.02 at fc = 0.1, 2.39 on 1.31 at 0.25).
//   - Flush (the plan's I6): window_frames() is a design-time constant, the
//     zero input frames after which the exact recursion's state is below the
//     smallest normal value of S from any state S can hold. With zero input
//     the cascade is x' = A x, A its 2(N/2) x 2(N/2) state-transition matrix;
//     the window is the smallest n with
//
//         ln(max_finite) + ln ||A^n||_inf  <  ln(min_normal),
//
//     plus a margin of one frame per section, found by binary lifting over
//     A^(2^j) kept normalized with a separate log scale (O(log n) matrix
//     products, even at a cutoff of a few hertz in double).
//     What the arithmetic leaves after it is the host's FTZ's
//     (iir_kernels.h's honest limit): usually exact silence; without FTZ a
//     subnormal residue; under FTZ, for some designs, a limit cycle just
//     above the smallest normal value (FTZ flushes a2 y while s1 survives,
//     and |a1| > 1 lifts it again), growing as the cutoff falls (measured on
//     x86-64, from full scale: up to 124 times std::numeric_limits<S>::min()
//     at order 8, fc = 0.002, float; 99 times at order 4, fc = 0.002, double;
//     none for most designs between 0.01 and 0.3). Pinned bound: 256 times
//     the smallest normal value (3e-36 in float, about -700 dB re full
//     scale), orders 4 and 8, either FTZ setting
//     (`butterworth_test.TheTailAfterTheWindowIsBounded` under FTZ,
//     `butterworth_subnormals.TheTailAfterTheWindowIsBounded` without, in
//     the tap_dsp_iir_subnormals executable). In a chain, flush() therefore
//     equals zero padding up to that residue rather than bit for bit, and
//     what follows it is the residue, not exact silence (chain.h).
//   - Non-finite input poisons the state until reset() (I9: no per-sample
//     check).
//   - Formats: double (golden) and float (float states and products, I14).
//     Float against double, order 4, white noise at 0.9 of full scale:
//     1.9e-7 at fc = 0.25, 4.3e-7 at fc = 0.07, 6.0e-6 at fc = 0.01 (x86-64;
//     transposed DF II's float error grows as the poles approach z = 1);
//     the battery bounds fc = 0.07 at 2e-6 (`TracksTheDoubleRecursion`).
//     Q15 / Q31: the kernel has them (direct form I, iir_kernels.h); this
//     stage does not offer them yet (a Butterworth at a low cutoff needs the
//     L-infinity section scaling iir_design.h does for the Chebyshev II).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "tap/dsp/detail/expects.h"
#include "tap/dsp/iir_design.h"
#include "tap/dsp/iir_kernels.h"

namespace tap::dsp {

    /// Whether a Butterworth cutoff is designable: 0 < fc < 1/2 cycles per sample, finite.
    constexpr bool butterworth_supports_cutoff(double fc) noexcept {
        return fc > 0.0 && fc < 0.5;
    }

    namespace detail {

        /// A square matrix kept as M * e^log_scale, M normalized to ||M||_inf = 1.
        template <std::size_t D>
        struct scaled_matrix {
            std::array<double, D * D> m{};
            double                    log_scale = 0.0;

            void normalize() noexcept {
                double n = 0.0;
                for (std::size_t r = 0; r < D; ++r) {
                    double row = 0.0;
                    for (std::size_t c = 0; c < D; ++c) {
                        row += std::fabs(m[r * D + c]);
                    }
                    n = std::max(n, row);
                }
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

        /// The zero-input state-transition matrix of a transposed direct form
        /// II cascade of K sections, state (s1, s2) per section in processing
        /// order.
        template <std::size_t K>
        scaled_matrix<2 * K> cascade_transition(const std::array<iir_section, K>& sec) noexcept {
            constexpr std::size_t D = 2 * K;
            scaled_matrix<D>      a;
            std::array<double, D> x{}; // section 0's input: the zero input, as a row over the old state
            for (std::size_t s = 0; s < K; ++s) {
                const iir_section&    c = sec[s];
                std::array<double, D> y = x; // y = b0 x + s1 (old)
                for (double& v : y) {
                    v *= c.b0;
                }
                y[2 * s] += 1.0;
                for (std::size_t j = 0; j < D; ++j) { // s1' = b1 x - a1 y + s2; s2' = b2 x - a2 y
                    a.m[(2 * s) * D + j]     = c.b1 * x[j] - c.a1 * y[j];
                    a.m[(2 * s + 1) * D + j] = c.b2 * x[j] - c.a2 * y[j];
                }
                a.m[(2 * s) * D + 2 * s + 1] += 1.0;
                x = y; // the next section's input is this one's output
            }
            a.normalize();
            return a;
        }

        /// The flush window (the file header's formula): zero input frames
        /// after which the exact recursion's state is below S's smallest
        /// normal value from any state S can hold, plus one frame per section.
        template <typename S, std::size_t K>
        std::size_t cascade_flush_window(const std::array<iir_section, K>& sections) {
            using mat          = scaled_matrix<2 * K>;
            const mat    a     = cascade_transition(sections);
            const double start = std::log(static_cast<double>(std::numeric_limits<S>::max()));
            const double goal  = std::log(static_cast<double>(std::numeric_limits<S>::min()));
            auto         done  = [&](const mat& p) { return start + p.log_scale < goal; };

            std::vector<mat> pow2{a}; // A^(2^j) until one alone is done
            while (!done(pow2.back()) && pow2.size() < 63) {
                pow2.push_back(pow2.back() * pow2.back());
            }
            mat         acc; // binary lifting: the largest n from those powers not yet done
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
            mat p = have ? acc * a : a;
            ++n;
            while (!done(p)) {
                p = p * a;
                ++n;
            }
            return n + K;
        }

    } // namespace detail

    /// The N/2 sections of an order-N Butterworth low-pass at fc cycles per
    /// sample, in increasing pole radius. Precondition: supports_cutoff(fc).
    template <std::size_t N>
    std::array<iir_section, N / 2> design_butterworth_lowpass(double fc) noexcept {
        static_assert(N >= 2 && N <= 16 && N % 2 == 0, "design_butterworth_lowpass: even orders 2 ... 16");
        TAP_EXPECTS(butterworth_supports_cutoff(fc));
        constexpr double               pi = std::numbers::pi;
        const double                   k  = std::tan(pi * fc);
        const double                   k2 = k * k;
        std::array<iir_section, N / 2> sec{};
        // k = N/2 - 1 has the largest alpha (lowest Q, smallest radius): run it first.
        for (std::size_t i = 0; i < N / 2; ++i) {
            const std::size_t idx = N / 2 - 1 - i;
            const double alpha = 2.0 * std::sin(static_cast<double>(2 * idx + 1) * pi / (2.0 * static_cast<double>(N)));
            const double a0    = 1.0 + alpha * k + k2;
            sec[i].b0          = k2 / a0;
            sec[i].b1          = 2.0 * k2 / a0;
            sec[i].b2          = k2 / a0;
            sec[i].a1          = 2.0 * (k2 - 1.0) / a0;
            sec[i].a2          = (1.0 - alpha * k + k2) / a0;
        }
        return sec;
    }

    /// An order-N Butterworth low-pass as a 1:1 chain stage over Ch
    /// interleaved channels, float or double (the file header's contract).
    template <typename S, std::size_t N, std::size_t Ch = 1>
    class basic_butterworth_lowpass {
        static_assert(N >= 2 && N <= 16 && N % 2 == 0, "basic_butterworth_lowpass: even orders 2 ... 16");
        static_assert(N / 2 <= k_max_cheby2_sections, "basic_butterworth_lowpass: the kernel's section capacity");
        static_assert(Ch >= 1, "basic_butterworth_lowpass: at least one channel");
        static_assert(std::is_same_v<S, float> || std::is_same_v<S, double>,
                      "basic_butterworth_lowpass: float and double (the file header: Q15 / Q31 not offered yet)");

      public:
        using sample                            = S;
        static constexpr std::size_t k_up       = 1;
        static constexpr std::size_t k_down     = 1;
        static constexpr std::size_t k_order    = N;
        static constexpr std::size_t k_sections = N / 2;
        static constexpr std::size_t channels() noexcept { return Ch; }

        /// The construction-time check, as a value: empty when cutoff is designable.
        static constexpr const char* validate(double cutoff) noexcept {
            return butterworth_supports_cutoff(cutoff)
                       ? nullptr
                       : "basic_butterworth_lowpass: cutoff must be in (0, 1/2) cycles per sample";
        }

        /// cutoff: the -3 dB frequency in cycles per sample, 0 < cutoff < 1/2.
        explicit basic_butterworth_lowpass(double cutoff)
            : m_cutoff(cutoff) {
            if (const char* why = validate(cutoff)) {
                raise<std::invalid_argument>(why);
            }
            m_sections           = design_butterworth_lowpass<N>(cutoff);
            m_table              = biquad_cascade_table<S>(std::span<const iir_section>(m_sections));
            m_window             = detail::cascade_flush_window<S>(m_sections);
            m_nominal            = sections_group_delay(m_sections, 0.0);
            double        dev    = 0.0;
            constexpr int k_grid = 1024;
            for (int i = 0; i <= k_grid; ++i) {
                const double f = cutoff * static_cast<double>(i) / k_grid;
                dev            = std::max(dev, std::fabs(sections_group_delay(m_sections, f) - m_nominal));
            }
            m_deviation = dev;
        }

        /// n interleaved frames in, n out; in and out may be the same buffer.
        /// Precondition (debug builds, where a bit governs S): the calling
        /// thread flushes subnormals (subnormals.h).
        std::size_t process(const S* in, std::size_t n, S* out) noexcept {
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t c = 0; c < Ch; ++c) {
                    out[i * Ch + c] = biquad_cascade<S>(m_table, m_state[c].data(), in[i * Ch + c]);
                }
            }
            return n;
        }

        std::size_t outputs_for(std::size_t n) const noexcept { return n; }
        void        reset() noexcept {
            for (auto& st : m_state) {
                st.fill(S(0));
            }
        }

        /// Zero input frames after which the exact recursion's state is below
        /// the smallest normal value of S (the residue after it: the file
        /// header's Flush).
        std::size_t window_frames() const noexcept { return m_window; }

        /// The group delay at f cycles per sample, 0 <= f < 1/2, in frames:
        /// -d(phase)/d(omega) of the double design.
        double group_delay_output_frames(double f) const noexcept { return sections_group_delay(m_sections, f); }

        /// The nominal latency: the group delay at DC, in frames. Delay a
        /// parallel path by this to align it with the stage's output.
        double nominal_latency_output_frames() const noexcept { return m_nominal; }

        /// The passband group-delay deviation, in frames: max over 0 <= f <=
        /// cutoff of |group_delay_output_frames(f) - nominal|, on a
        /// 1,024-point grid. The alignment error left after compensating the
        /// nominal alone (chain.h's reports_nominal_latency); 0 would mean the
        /// passband is delayed alike. A Butterworth's peaks at the band edge.
        double latency_deviation_output_frames() const noexcept { return m_deviation; }

        double                                cutoff() const noexcept { return m_cutoff; }
        const std::array<iir_section, N / 2>& sections() const noexcept { return m_sections; }

        /// True when every channel's state is exactly zero (silence from here on).
        bool is_silent() const noexcept {
            for (const auto& st : m_state) {
                for (const S v : st) {
                    if (v != S(0)) {
                        return false;
                    }
                }
            }
            return true;
        }

      private:
        double                                      m_cutoff;
        std::array<iir_section, N / 2>              m_sections{};
        biquad_cascade_table<S>                     m_table{};
        std::array<std::array<iir_state<S>, N>, Ch> m_state{}; ///< z1 then z2 per section (iir_kernels.h)
        std::size_t                                 m_window    = 0;
        double                                      m_nominal   = 0.0;
        double                                      m_deviation = 0.0;
    };

    template <std::size_t Ch = 1>
    using butterworth4 = basic_butterworth_lowpass<double, 4, Ch>; ///< order 4, the golden model
    template <std::size_t Ch = 1>
    using butterworth4_32 = basic_butterworth_lowpass<float, 4, Ch>; ///< order 4, the embedded profile

} // namespace tap::dsp
