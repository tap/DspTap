/// @file decimate.h
/// @brief Fixed integer-ratio decimators to 16 kHz (by 2, 3, 6), on the FIR substrate.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The host-rate stage in front of a 16 kHz speech front end: 32, 48 or
// 96 kHz in, 16 kHz out. Built in RatioTap's pattern — the ratio is a
// compile-time type, the prototype is a Kaiser-windowed sinc from kaiser.h,
// the hot loop is fir_kernels.h's dot_row over the sample_traits.h formats
// (float golden, Q15/Q31 fixed point) — but deliberately NOT RatioTap: that
// library's charter is 44.1 <-> 48 only. A 44.1 kHz host composes RatioTap's
// 44.1 -> 48 in front of the by-3 stage here.
//
// Contract, as numbers:
//   - Ratios: 2, 3 and 6 (input 32 / 48 / 96 kHz for a 16 kHz output). The
//     static_assert is the charter; widening it is a documented decision.
//   - Output alignment: y[k] = sum_t h[t] x[k*M - t], with x[n] = 0 for n < 0.
//     The k-th output is produced as input sample k*M arrives, so a stream of
//     n inputs yields ceil(n / M) outputs from a fresh instance; outputs_for()
//     gives the exact count from the current phase.
//   - Prototype: Kaiser-windowed sinc, cutoff at the OUTPUT Nyquist
//     (cutoff_norm = 1/M in kaiser.h's input-Nyquist units), sum(h) = 1
//     (DC gain 1), odd tap count so the group delay is the integer
//     (taps - 1) / 2 input samples; linear phase.
//   - Profiles (stopband attenuation, flat passband, taps by 2 / 3 / 6) — the
//     tap counts are the minimal odd counts meeting the stopband with >= 1 dB
//     margin and <= 0.1 dB passband deviation on a 25 Hz grid, searched by
//     tools/reference/make_frontend_reference.py and pinned by the tests:
//
//       | profile     | stopband | passband | by 2 | by 3 | by 6 |
//       |-------------|----------|----------|------|------|------|
//       | economy     |  70 dB   | 7000 Hz  |  81  | 121  | 239  |
//       | transparent | 100 dB   | 7600 Hz  | 259  | 389  | 773  |
//
//     The stopband edge is 16000 - passband, so nothing above the passband
//     edge can alias into it. economy is the default: a keyword spotter's mel
//     bands stop at 7.6 kHz and its features carry little above 7 kHz, and
//     economy costs 121 MACs per 16 kHz output (by 3) — under 20k MACs per
//     10 ms hop. transparent exists for offline use.
//   - Sample formats: float (double accumulation, the golden model, pinned
//     sample-for-sample against a committed numpy reference), Q15 and Q31
//     through sample_traits.h with row-sum-preserving quantization so DC
//     gain stays exactly 1. Mono: the consumer is single-channel by charter.
//
// Construction designs the filter (runtime double, off the audio path) and
// allocates; process() and reset() are noexcept and allocation-free.

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <span>
#include <vector>

#include "tap/dsp/fir_kernels.h"
#include "tap/dsp/kaiser.h"
#include "tap/dsp/quantize.h"
#include "tap/dsp/sample_traits.h"

namespace tap::dsp {

    /// Quality profile for the decimators; tap counts pinned per ratio.
    struct decimate_profile {
        double      passband_hz       = 7000.0;
        double      stopband_atten_db = 70.0;
        std::size_t taps_by_2         = 81;
        std::size_t taps_by_3         = 121;
        std::size_t taps_by_6         = 239;

        /// The default: 70 dB, flat to 7 kHz.
        static constexpr decimate_profile economy() noexcept { return {}; }

        /// 100 dB, flat to 7.6 kHz — the offline tier.
        static constexpr decimate_profile transparent() noexcept {
            return {.passband_hz       = 7600.0,
                    .stopband_atten_db = 100.0,
                    .taps_by_2         = 259,
                    .taps_by_3         = 389,
                    .taps_by_6         = 773};
        }

        template <std::size_t M>
        constexpr std::size_t taps() const noexcept {
            if constexpr (M == 2) {
                return taps_by_2;
            }
            else if constexpr (M == 3) {
                return taps_by_3;
            }
            else {
                return taps_by_6;
            }
        }
    };

    /// Compile-time facts of one ratio.
    template <std::size_t M>
    struct decimate_traits {
        static_assert(M == 2 || M == 3 || M == 6, "decimate.h serves the 16 kHz front end: ratios 2, 3 and 6 only");
        static constexpr std::size_t k_ratio          = M;
        static constexpr double      k_output_rate_hz = 16000.0;
        static constexpr double      k_input_rate_hz  = 16000.0 * static_cast<double>(M);
    };

    /// Designs the profile's prototype for ratio M: `taps` coefficients at the
    /// input rate, sum 1, cutoff at the output Nyquist. Allocates.
    template <std::size_t M>
    inline std::vector<double> design_decimator(const decimate_profile& p) {
        std::vector<double> h(p.taps<M>());
        design_prototype(h, 1, 1.0 / static_cast<double>(M), kaiser_beta(p.stopband_atten_db));
        return h;
    }

    /// Fixed-ratio decimator over one sample format.
    template <sample_type S, std::size_t M>
    class basic_decimator {
      public:
        using traits = decimate_traits<M>;
        using coeff  = typename sample_traits<S>::coeff;

        static constexpr std::size_t k_ratio = M;

        explicit basic_decimator(const decimate_profile& p = decimate_profile::economy())
            : m_taps(p.taps<M>()) {
            assert(m_taps >= 3 && (m_taps % 2) == 1);
            const std::vector<double> proto = design_decimator<M>(p);
            m_h.resize(m_taps);
            quantize_row_preserving_sum<S>(proto, m_h);
            m_buf.assign(m_taps - 1 + k_block, sample_traits<S>::silence());
            m_fill = m_taps - 1;
        }

        /// Decimate in_frames samples; writes outputs_for(in_frames) samples.
        /// Bit-identical for any chunking of the same stream.
        std::size_t process(const S* in, std::size_t in_frames, S* out) noexcept {
            std::size_t made = 0;
            for (std::size_t i = 0; i < in_frames; ++i) {
                if (m_fill == m_buf.size()) {
                    std::memmove(m_buf.data(), m_buf.data() + (m_fill - (m_taps - 1)), (m_taps - 1) * sizeof(S));
                    m_fill = m_taps - 1;
                }
                m_buf[m_fill++] = in[i];
                if (m_phase == 0) {
                    // Newest sample is x[k*M]; the window holds x[k*M - taps + 1 .. k*M],
                    // oldest first, and the symmetric prototype makes the forward dot
                    // equal to the convolution sum exactly (h[t] == h[taps - 1 - t]).
                    out[made++] = dot_row<S>(m_h.data(), m_buf.data() + (m_fill - m_taps), m_taps);
                }
                m_phase = (m_phase + 1) % M;
            }
            return made;
        }

        /// Outputs process() will write for in_frames more inputs.
        std::size_t outputs_for(std::size_t in_frames) const noexcept {
            // The j-th next input (j from 1) emits when (m_phase + j - 1) % M == 0,
            // so the first emitting input is j = 1 + (M - m_phase) % M.
            const std::size_t first = 1 + (M - m_phase) % M;
            if (in_frames < first) {
                return 0;
            }
            return 1 + (in_frames - first) / M;
        }

        void reset() noexcept {
            std::fill(m_buf.begin(), m_buf.end(), sample_traits<S>::silence());
            m_fill  = m_taps - 1;
            m_phase = 0;
        }

        std::size_t            taps() const noexcept { return m_taps; }
        std::span<const coeff> coefficients() const noexcept { return m_h; }
        /// Group delay in input samples: (taps - 1) / 2, an integer.
        std::size_t latency_input_samples() const noexcept { return (m_taps - 1) / 2; }

      private:
        static constexpr std::size_t k_block = 256; ///< history buffer slack before the memmove

        std::size_t        m_taps;
        std::vector<coeff> m_h;
        std::vector<S>     m_buf; ///< taps - 1 history samples followed by up to k_block new ones
        std::size_t        m_fill  = 0;
        std::size_t        m_phase = 0; ///< inputs since the last emitted output, in [0, M)
    };

    using decimate_by_2 = basic_decimator<float, 2>; ///< 32 kHz -> 16 kHz, float golden
    using decimate_by_3 = basic_decimator<float, 3>; ///< 48 kHz -> 16 kHz
    using decimate_by_6 = basic_decimator<float, 6>; ///< 96 kHz -> 16 kHz

} // namespace tap::dsp
