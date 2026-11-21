/// @file butterworth.h
/// @brief Butterworth low-pass: the even-order design (bilinear, pre-warped) and a 1:1 chain stage over it.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The first slice of SampleRateTap's IIR plan's DspTap part (IIR-a): the
// biquad cascade kernel (biquad.h), the chain's nominal latency (chain.h) and
// one recursive stage that exercises them, a Butterworth low-pass, under the
// plan's contracts for a recursive stage in a chain: latency I5 and flush
// I6. Subnormals are the client's (I9: the client sets flush-to-zero; this
// code neither sets it nor checks it).
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
//     (supports_cutoff; the constructor rejects anything else through
//     raise<std::invalid_argument>). Magnitude |H(fc)| = 1/sqrt(2)
//     (-3.0103 dB), |H(0)| = 1, H(1/2) = 0; monotone in between.
//   - A 1:1 sync_stage (chain.h): k_up = k_down = 1, channels() = Ch, any
//     Ch >= 1, frames interleaved; process() writes exactly n frames
//     (outputs_for(n) = n), noexcept, allocation-free, bit-identical for any
//     chunking (the state is per channel and per sample); across channel
//     counts and chain.h's flush and process paths, given fft.h's
//     fp-contraction policy (contraction off; fft_fir.h says why).
//   - Latency (the plan's I5): a recursive stage has no single delay, so it
//     never provides latency_output_frames(). It reports
//     group_delay_output_frames(f) (biquad.h's formula, double, at f cycles
//     per sample, 0 <= f < 1/2), nominal_latency_output_frames() = the
//     group delay at DC, and passband_group_delay_deviation_output_frames()
//     = the largest |group_delay(f) - nominal| over the passband 0 <= f <=
//     fc, on a 1,024-point grid (a Butterworth's group delay peaks near fc).
//   - Flush (the plan's I6): window_frames() = biquad.h's flush window, a
//     design-time constant: the zero input frames after which the exact
//     recursion's state is below the smallest normal value. What the
//     arithmetic leaves after it depends on the client (biquad.h): usually
//     exact silence; without FTZ a subnormal residue (below the smallest
//     normal value); under FTZ, for some designs, a limit cycle just above
//     it, growing as the cutoff falls (measured on x86-64, from full scale:
//     up to 124 times std::numeric_limits<S>::min() at order 8, fc = 0.002,
//     float; 99 times at order 4, fc = 0.002, double; none for most designs
//     between 0.01 and 0.3). Pinned bound: 256 times the smallest normal
//     value (3e-36 in float, about -700 dB re full scale), orders 2 ... 8,
//     cutoffs 0.002 ... 0.49, either FTZ setting
//     (`butterworth_test.TheTailAfterTheWindowIsBounded`). In a chain, flush()
//     therefore equals zero padding up to that residue rather than bit for
//     bit, and what follows it is the residue, not exact silence (chain.h).
//   - Non-finite input poisons the state until reset() (I9: no per-sample
//     check).
//   - Formats: double (golden) and float (float states and products, I14).
//     Float against double, order 4, white noise at 0.9 of full scale:
//     1.9e-7 at fc = 0.25, 4.3e-7 at fc = 0.07, 6.0e-6 at fc = 0.01 (x86-64;
//     transposed DF II's float error grows as the poles approach z = 1);
//     the battery bounds fc = 0.07 at 2e-6 (`TracksTheDoubleRecursion`).
//     Q15 / Q31 wait for the plan's I8 (biquad.h).

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "tap/dsp/biquad.h"
#include "tap/dsp/detail/expects.h"

namespace tap::dsp {

    /// Whether a Butterworth cutoff is designable: 0 < fc < 1/2 cycles per sample, finite.
    constexpr bool butterworth_supports_cutoff(double fc) noexcept {
        return fc > 0.0 && fc < 0.5;
    }

    /// The N/2 sections of an order-N Butterworth low-pass at fc cycles per
    /// sample, in increasing pole radius. Precondition: supports_cutoff(fc).
    template <std::size_t N>
    std::array<biquad_coefficients, N / 2> design_butterworth_lowpass(double fc) noexcept {
        static_assert(N >= 2 && N <= 16 && N % 2 == 0, "design_butterworth_lowpass: even orders 2 ... 16");
        TAP_EXPECTS(butterworth_supports_cutoff(fc));
        constexpr double                       pi = 3.14159265358979323846;
        const double                           k  = std::tan(pi * fc);
        const double                           k2 = k * k;
        std::array<biquad_coefficients, N / 2> sec{};
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

    /// An order-N Butterworth low-pass as a 1:1 chain stage over Ch interleaved channels.
    template <typename S, std::size_t N, std::size_t Ch = 1>
    class basic_butterworth_lowpass {
        static_assert(N >= 2 && N <= 16 && N % 2 == 0, "basic_butterworth_lowpass: even orders 2 ... 16");
        static_assert(Ch >= 1, "basic_butterworth_lowpass: at least one channel");

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
            m_sections = design_butterworth_lowpass<N>(cutoff);
            m_kernel.set(m_sections);
            m_window             = flush_window_frames<S>(m_sections);
            m_nominal            = group_delay_samples(m_sections, 0.0);
            double        dev    = 0.0;
            constexpr int k_grid = 1024;
            for (int i = 0; i <= k_grid; ++i) {
                const double f = cutoff * static_cast<double>(i) / k_grid;
                dev            = std::max(dev, std::fabs(group_delay_samples(m_sections, f) - m_nominal));
            }
            m_deviation = dev;
        }

        /// n interleaved frames in, n out.
        std::size_t process(const S* in, std::size_t n, S* out) noexcept {
            m_kernel.process(in, n, out);
            return n;
        }

        std::size_t outputs_for(std::size_t n) const noexcept { return n; }
        void        reset() noexcept { m_kernel.reset(); }

        /// Zero input frames after which the exact recursion's state is below
        /// the smallest normal value (the residue: the file header).
        std::size_t window_frames() const noexcept { return m_window; }

        /// Group delay at f cycles per sample, 0 <= f < 1/2, in frames.
        double group_delay_output_frames(double f) const noexcept { return group_delay_samples(m_sections, f); }
        /// The group delay at DC, in frames.
        double nominal_latency_output_frames() const noexcept { return m_nominal; }
        /// max |group delay - nominal| over 0 <= f <= cutoff, in frames.
        double passband_group_delay_deviation_output_frames() const noexcept { return m_deviation; }

        double                                        cutoff() const noexcept { return m_cutoff; }
        const std::array<biquad_coefficients, N / 2>& sections() const noexcept { return m_sections; }
        /// True when the state is exactly zero.
        bool is_silent() const noexcept { return m_kernel.is_silent(); }

      private:
        double                                 m_cutoff;
        std::array<biquad_coefficients, N / 2> m_sections{};
        basic_biquad_cascade<S, N / 2, Ch>     m_kernel;
        std::size_t                            m_window    = 0;
        double                                 m_nominal   = 0.0;
        double                                 m_deviation = 0.0;
    };

    template <std::size_t Ch = 1>
    using butterworth4 = basic_butterworth_lowpass<double, 4, Ch>; ///< order 4, the golden model
    template <std::size_t Ch = 1>
    using butterworth4_32 = basic_butterworth_lowpass<float, 4, Ch>; ///< order 4, the embedded profile

} // namespace tap::dsp
