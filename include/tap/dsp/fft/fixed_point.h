/// @file fixed_point.h
/// @brief The fixed-point real FFT engine: one int32 radix-4 kernel for the Q15 and Q31 profiles.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Stage 3b of docs/audit-fft-and-code-smells.md (design: Part 7; the
// arithmetic it is written against: fft/fft_arith.h, whose docstrings are
// this kernel's specification). basic_real_fft<std::int16_t, Scaling> and
// basic_real_fft<std::int32_t, Scaling> in fft.h route here; the class has
// the same engine shape as the floating engine (constructor from n,
// forward_inplace, inverse_inplace, size) so Stage 4 can make engines
// explicit.
//
// The contract this engine shares with the floating profiles is the packing
// (a[0] DC, a[1] Nyquist, a[2k] + i a[2k+1] for 1 <= k < N/2), the sign
// convention W = exp(+2*pi*i/N), and the unnormalized inverse; what it adds
// is an exponent. Every transform returns e such that the double golden
// model's result on the same input, read as fractions of full scale, equals
// this buffer times 2^e (up to rounding noise). Details, per profile and
// policy, in the class docstring below and in docs/fft-fixed-point.md.
//
// Implemented from the published literature only (house IP policy):
//   P. D. Welch, "A fixed-point fast Fourier transform error analysis,"
//     IEEE Trans. Audio Electroacoust. AU-17(2), 1969 — the per-stage
//     scaling analysis and the noise model the floors are checked against.
//   A. V. Oppenheim and C. J. Weinstein, "Effects of finite register length
//     in digital filtering and the fast Fourier transform," Proc. IEEE
//     60(8), 1972 — rounding-noise accumulation across stages.
//   A. V. Oppenheim and R. W. Schafer, Discrete-Time Signal Processing, 3rd
//     ed., Sec. 9.7 — block floating point.
//   The real post-pass formulas and the DC/Nyquist glue are transcribed from
//   Ooura's rdft/rftfsub/rftbsub (third_party/ooura/fftsg.c) into the
//   trait's arithmetic; the complex kernel is not Ooura's (Decision D2).

#pragma once

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

#include "tap/dsp/fft/fft_arith.h"
#include "tap/dsp/fft/tables.h"

namespace tap::dsp {

    /// Scaling policies of the fixed-point real FFT profiles (Decision D3).
    /// A template parameter of basic_real_fft for std::int16_t / std::int32_t
    /// samples; the floating profiles accept only scaling::fixed.
    namespace scaling {
        /// Shift-before-butterfly by the stage's growth every time: 2 bits per
        /// radix-4 stage, 1 per radix-2, 1 before the real post-pass, plus the
        /// profile's input pre-shift. The exponent is a constant of N
        /// (fixed_scaling_exponent), the forward output is exactly X / N in
        /// the sample's Q format for Q15 and X / 2N for Q31, and the noise
        /// floor is the same for every input. For magnitude and power
        /// consumers (log-mel front ends, detectors).
        struct fixed {};

        /// Before each stage, a headroom scan of the block; the stage shifts
        /// only as much as its growth requires beyond the headroom present
        /// (never consuming the last bit), and the exponent accumulates what
        /// was shifted: 0 <= e <= fixed_scaling_exponent(N), data-dependent.
        /// A quiet block keeps its precision instead of losing log2 N bits.
        /// For round-trip consumers (overlap-add resynthesis, adaptive filters).
        struct block_floating {};
    } // namespace scaling

    namespace detail {

        /// The fixed-point real FFT engine over the fft_arith<Sample> trait.
        ///
        /// Sample is std::int16_t (Q0.15 I/O, computed in int32 as Q2.29 after
        /// fft_arith<std::int16_t>::widen) or std::int32_t (Q0.31, computed in
        /// place). Scaling is scaling::fixed or scaling::block_floating.
        ///
        /// Structure (audit Part 7): the N real samples are read as N/2
        /// complex values z[j] = x[2j] + i x[2j+1]; a radix-4
        /// decimation-in-frequency kernel of length M = N/2 (radix-2 final
        /// stage when log2 M is odd) with outputs placed so the result is
        /// bit-reversed, the permutation from make_bit_reversal_table, then
        /// Ooura's real post-pass (rftfsub) and DC/Nyquist glue transcribed
        /// into the trait's operations. The inverse mirrors rdft's isgn < 0
        /// path: the DC/Nyquist pre-glue, rftbsub, the conjugate kernel.
        /// Every multiply is fft_arith<>::work::mul_coeff by a Q1.30 twiddle
        /// (one rounding per real product, two per complex product), every
        /// sum is the saturating add/sub, every scale is shr_round; the kernel
        /// invents no arithmetic of its own.
        ///
        /// Exponent contract. Let G be basic_real_fft<double> on the same
        /// input read as fractions of full scale. forward_inplace leaves the
        /// packed spectrum such that G.forward_inplace's result == data * 2^e
        /// and returns e; inverse_inplace leaves G's UNNORMALIZED
        /// inverse_inplace result == data * 2^e and returns e. No 2/N is ever
        /// applied: the exponent carries the scale, and a round trip gives
        /// x == out * 2^(e_fwd + e_inv + 1 - log2 N) under both policies
        /// (Ooura's unnormalized round trip has gain N/2).
        ///   - scaling::fixed: e == fixed_scaling_exponent(N)
        ///     == log2 N + fft_arith<Sample>::k_fixed_scaling_input_pre_shift
        ///     in both directions (Q15: log2 N, the forward output is exactly
        ///     X / N; Q31: log2 N + 1, exactly X / 2N).
        ///   - scaling::block_floating: 0 <= e <= fixed_scaling_exponent(N),
        ///     data-dependent; before each stage the shift is the smaller of
        ///     what the headroom scan demands (growth + 1 - headroom, so one
        ///     bit is never consumed) and what the fixed policy would have
        ///     accumulated by that stage, which is what makes the upper bound
        ///     a contract rather than a typical value. An all-zero block
        ///     returns e == 0.
        ///
        /// Size: 4 <= N <= 65536, a power of two, fixed at construction
        /// (asserted, as for the floating profiles). Everything is allocated
        /// in the constructor: the three tables, and for Q15 an int32 work
        /// buffer of N (the in-place API is preserved at the caller's int16
        /// buffer; Q31 transforms in place with no work buffer). The
        /// transforms are noexcept and allocation-free; the object is
        /// copyable; there is no alignment requirement; one transform at a
        /// time per object.
        template <typename Sample, typename Scaling>
        class fixed_point_rdft {
            static_assert(std::is_same_v<Sample, std::int16_t> || std::is_same_v<Sample, std::int32_t>,
                          "fixed_point_rdft: the fixed-point profiles are std::int16_t (Q15) and std::int32_t (Q31)");
            static_assert(std::is_same_v<Scaling, scaling::fixed> || std::is_same_v<Scaling, scaling::block_floating>,
                          "fixed_point_rdft: Scaling is scaling::fixed or scaling::block_floating");

          public:
            using sample = Sample;
            using arith  = fft_arith<Sample>;
            using wide   = typename arith::wide;  ///< int32: the kernel's data width
            using work   = typename arith::work;  ///< the int32 butterfly arithmetic
            using coeff  = typename arith::coeff; ///< Q1.30 twiddle

            static constexpr bool        k_block_floating = std::is_same_v<Scaling, scaling::block_floating>;
            static constexpr std::size_t k_min_size       = 4;
            static constexpr std::size_t k_max_size       = 65536;

            /// The constant exponent of scaling::fixed, and the upper bound of
            /// scaling::block_floating: log2 n + the profile's input pre-shift
            /// (0 for Q15, whose two guard bits cover the growth; 1 for Q31).
            /// @pre n is a power of two in [k_min_size, k_max_size].
            [[nodiscard]] static constexpr int fixed_scaling_exponent(std::size_t n) noexcept {
                return log2_of(n) + arith::k_fixed_scaling_input_pre_shift;
            }

            /// Builds the tables (and the Q15 work buffer) for size n.
            /// @pre n is a power of two in [k_min_size, k_max_size].
            explicit fixed_point_rdft(std::size_t n)
                : m_n(n)
                , m_log2_n(log2_of(n))
                , m_twiddles(make_twiddle_table(n / 2))
                , m_post(make_real_post_pass_table(n))
                , m_bitrev(make_bit_reversal_table(n / 2))
                , m_work(std::is_same_v<Sample, std::int16_t> ? n : 0, wide{0}) {
                assert(n >= k_min_size && n <= k_max_size && (n & (n - 1)) == 0);
            }

            [[nodiscard]] std::size_t size() const noexcept { return m_n; }

            /// In-place forward transform: N samples -> packed spectrum, scaled
            /// by 2^-e. @return e (see the class docstring).
            int forward_inplace(Sample* data) noexcept {
                // TODO(3b): kernel. Skeleton: data unchanged, fixed exponent.
                (void)data;
                return fixed_scaling_exponent(m_n);
            }

            /// In-place inverse transform: packed spectrum -> N samples, the
            /// unnormalized inverse scaled by 2^-e. @return e.
            int inverse_inplace(Sample* data) noexcept {
                // TODO(3b): kernel. Skeleton: data unchanged, fixed exponent.
                (void)data;
                return fixed_scaling_exponent(m_n);
            }

          private:
            static constexpr int log2_of(std::size_t n) noexcept { return std::bit_width(n) - 1; }

            std::size_t                m_n;
            int                        m_log2_n;
            std::vector<coeff>         m_twiddles; ///< W_M^k, M = N/2, interleaved (tables.h)
            std::vector<coeff>         m_post;     ///< real post-pass (wkr, wki) per bin
            std::vector<std::uint32_t> m_bitrev;   ///< permutation of the M complex outputs
            std::vector<wide>          m_work;     ///< Q15 only: the widened block
        };

    } // namespace detail

} // namespace tap::dsp
