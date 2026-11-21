/// @file fft_fir.h
/// @brief An FIR filter by fast convolution (overlap-save) as a 1:1 chain stage over the real FFT.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Overlap-save (Oppenheim & Schafer, "Discrete-Time Signal Processing",
// 3rd ed., 2010, sec. 8.7.3; Stockham, "High-speed convolution and
// correlation", AFIPS 1966): the input is cut into blocks of B frames; each
// block and the one before it (N = 2B samples per channel) go through the
// real FFT, are multiplied by H, the spectrum of the L taps zero-padded to
// N, and come back through the inverse; the last B samples of the circular
// result equal the linear convolution, because B >= L - 1.
//
// Contract, as numbers:
//   - Taps: L >= 1, given in double at construction. Block B: a power of two
//     with B >= L - 1 and the FFT engine supporting N = 2B
//     (basic_real_fft<S>::supports_size); default B = the smallest power of
//     two >= max(L - 1, 16). validate(L, B) returns the reason the
//     constructor would reject through raise<std::invalid_argument>.
//   - A 1:1 sync_stage (chain.h): k_up = k_down = 1, channels() = Ch, any Ch
//     >= 1, frames interleaved; process() writes exactly n frames
//     (outputs_for(n) = n), noexcept, allocation-free, bit-identical for any
//     chunking (the blocks fall on the stream's own frame count, not the
//     call's). Bit-identical as well across channel counts (a channel of a
//     Ch-channel stage against a mono one) and across chain.h's flush and
//     process paths, given fft.h's fp-contraction policy (contraction off):
//     g++ fuses by default where the target has FMA, and may fuse two
//     inlined copies of one loop differently.
//   - Output: y_out[t] = sum_k h[k] x[t - B - k], with x[t] = 0 for t < 0:
//     the direct convolution delayed by exactly block_delay_frames() = B
//     frames (the first B outputs are silence). Accuracy against the direct
//     form in double: the FFT's (fft.h); measured 1.2e-7 of full scale in
//     float and 6e-16 in double for a 255-tap low-pass over white noise at
//     0.9 (x86-64, srdif), bounded at 1e-6 and 1e-14 by the battery on every
//     engine (srdif, CMSIS-DSP on the M55, vDSP on macOS)
//     (`fft_fir_test.MatchesDirectConvolution`).
//   - Latency: latency_output_frames() = B + (L - 1) / 2, exact, for taps
//     that are symmetric or antisymmetric (linear phase; linear_phase(), the
//     taps compared exactly at construction): the group delay. For other
//     taps the stage has no single delay and the call is a precondition
//     violation (TAP_EXPECTS(linear_phase())); block_delay_frames() is
//     always exact.
//   - Flush: window_frames() = 3 B, a constant. Zeros reach every sample of
//     an FFT window two whole blocks after the block they start in, and a
//     window of zeros transforms to exactly zero, so after 3 B zero frames
//     the output is exactly silence (chain.h's flush contract unchanged).
//   - Work: all of it happens on the frame that completes a block: 2 Ch
//     transforms of N plus N/2 + 1 complex products per channel; the other
//     B - 1 frames only copy. So the worst single call costs a block's
//     work, and a host's per-call budget is that, not the average.
//   - Layout: the FFT wants each channel contiguous, so a block
//     de-interleaves its input and re-interleaves its output. On x86-64 that
//     still costs less than one stage per channel (planar 1.04 to 1.06
//     times the instructions at 2 to 12 channels); the per-core counts are
//     bench/README.md's "Layout".
//   - Formats: double and float (the FFT's floating profiles; H in the
//     stage's format, computed in double and rounded). Q15 / Q31 are not
//     offered: the fixed-point FFT returns an exponent per transform, which
//     a streaming stage would have to carry across blocks.
//   - Lives in the FFT's ABI-tag namespace (fft.h, "The ABI tag"): it holds
//     a basic_real_fft<S> by value.

#pragma once

#include <algorithm>
#include <bit>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "tap/dsp/chain.h"
#include "tap/dsp/detail/expects.h"
#include "tap/dsp/fft.h"

namespace tap::dsp::inline TAP_DSP_FFT_ABI {

    /// An FIR filter by overlap-save fast convolution, Ch interleaved channels.
    template <typename S, std::size_t Ch = 1>
    class basic_fft_fir {
        static_assert(std::is_same_v<S, float> || std::is_same_v<S, double>,
                      "basic_fft_fir: float and double (the FFT's floating profiles)");
        static_assert(Ch >= 1, "basic_fft_fir: at least one channel");

      public:
        using sample                        = S;
        static constexpr std::size_t k_up   = 1;
        static constexpr std::size_t k_down = 1;
        static constexpr std::size_t channels() noexcept { return Ch; }

        /// The default block for L taps: the smallest power of two >= max(L - 1, 16).
        static constexpr std::size_t default_block(std::size_t taps) noexcept {
            return std::bit_ceil(std::max<std::size_t>(taps > 0 ? taps - 1 : 0, 16));
        }

        /// The construction-time check, as a value: nullptr when (taps, block) is accepted.
        static constexpr const char* validate(std::size_t taps, std::size_t block) noexcept {
            if (taps == 0) {
                return "basic_fft_fir: at least one tap";
            }
            if (block == 0 || !std::has_single_bit(block)) {
                return "basic_fft_fir: the block must be a power of two";
            }
            if (block < taps - 1) {
                return "basic_fft_fir: the block must be at least taps - 1";
            }
            if (!basic_real_fft<S>::supports_size(2 * block)) {
                return "basic_fft_fir: the FFT engine does not support twice the block";
            }
            return nullptr;
        }

        /// taps: the impulse response; block: 0 for default_block(taps.size()).
        explicit basic_fft_fir(std::span<const double> taps, std::size_t block = 0)
            : m_taps(taps.size())
            , m_block(block != 0 ? block : default_block(taps.size()))
            , m_fft(validate(m_taps, m_block) == nullptr ? 2 * m_block : basic_real_fft<S>::k_min_size) {
            if (const char* why = validate(m_taps, m_block)) {
                raise<std::invalid_argument>(why);
            }
            const std::size_t n = 2 * m_block;
            // H in double through the golden FFT, the 2/N of the inverse folded in.
            {
                basic_real_fft<double> fft64(n);
                std::vector<double>    h(n, 0.0);
                std::copy(taps.begin(), taps.end(), h.begin());
                fft64.forward_inplace(h.data());
                const double scale = 2.0 / static_cast<double>(n);
                m_h.resize(n);
                for (std::size_t i = 0; i < n; ++i) {
                    m_h[i] = static_cast<S>(h[i] * scale);
                }
            }
            bool sym = true, anti = true;
            for (std::size_t k = 0; k < m_taps; ++k) {
                sym  = sym && taps[k] == taps[m_taps - 1 - k];
                anti = anti && taps[k] == -taps[m_taps - 1 - k];
            }
            m_linear_phase = sym || anti;
            m_hist.assign(Ch * n, S(0));
            m_out.assign(Ch * m_block, S(0));
            m_work.assign(n, S(0));
        }

        /// n interleaved frames in, n out (delayed by block_delay_frames()).
        /// in and out may be the same buffer.
        std::size_t process(const S* in, std::size_t n, S* out) noexcept {
            const std::size_t b = m_block;
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t ch = 0; ch < Ch; ++ch) {
                    const S x                      = in[i * Ch + ch];
                    out[i * Ch + ch]               = m_out[m_pos * Ch + ch];
                    m_hist[ch * 2 * b + b + m_pos] = x;
                }
                if (++m_pos == b) {
                    run_block();
                    m_pos = 0;
                }
            }
            return n;
        }

        std::size_t outputs_for(std::size_t n) const noexcept { return n; }

        void reset() noexcept {
            std::fill(m_hist.begin(), m_hist.end(), S(0));
            std::fill(m_out.begin(), m_out.end(), S(0));
            m_pos = 0;
        }

        /// Zero input frames after which the output is exactly silence.
        std::size_t window_frames() const noexcept { return 3 * m_block; }

        /// The frames every output lags the direct convolution: B.
        std::size_t block_delay_frames() const noexcept { return m_block; }

        /// B + (L - 1) / 2, the group delay. Precondition: linear_phase().
        exact_ratio latency_output_frames() const noexcept {
            TAP_EXPECTS(m_linear_phase);
            return exact_ratio{2 * m_block + m_taps - 1, 2};
        }

        /// Whether the taps are symmetric or antisymmetric (exactly).
        bool        linear_phase() const noexcept { return m_linear_phase; }
        std::size_t taps() const noexcept { return m_taps; }
        std::size_t block() const noexcept { return m_block; }
        std::size_t fft_size() const noexcept { return 2 * m_block; }

      private:
        /// One block per channel: transform the 2B window, multiply by H,
        /// transform back, keep the last B; slide the window by B.
        void run_block() noexcept {
            const std::size_t b = m_block;
            const std::size_t n = 2 * b;
            for (std::size_t ch = 0; ch < Ch; ++ch) {
                S* hist = m_hist.data() + ch * n;
                std::copy(hist, hist + n, m_work.begin());
                m_fft.forward_inplace(m_work.data());
                S*       w = m_work.data();
                const S* h = m_h.data();
                w[0] *= h[0]; // DC
                w[1] *= h[1]; // Nyquist
                for (std::size_t k = 2; k < n; k += 2) {
                    const S re = w[k] * h[k] - w[k + 1] * h[k + 1];
                    const S im = w[k] * h[k + 1] + w[k + 1] * h[k];
                    w[k]       = re;
                    w[k + 1]   = im;
                }
                m_fft.inverse_inplace(w);
                for (std::size_t k = 0; k < b; ++k) {
                    m_out[k * Ch + ch] = w[b + k];
                }
                std::copy(hist + b, hist + n, hist);
            }
        }

        std::size_t       m_taps;
        std::size_t       m_block;
        basic_real_fft<S> m_fft;
        std::vector<S>    m_h;    ///< packed spectrum of the taps, times 2/N
        std::vector<S>    m_hist; ///< per channel: the previous block, then the one filling
        std::vector<S>    m_out;  ///< the block being emitted, interleaved
        std::vector<S>    m_work; ///< one channel's transform
        std::size_t       m_pos          = 0;
        bool              m_linear_phase = false;
    };

    template <std::size_t Ch = 1>
    using fft_fir = basic_fft_fir<double, Ch>; ///< the golden model
    template <std::size_t Ch = 1>
    using fft_fir32 = basic_fft_fir<float, Ch>; ///< the embedded profile

} // namespace tap::dsp::inline TAP_DSP_FFT_ABI
