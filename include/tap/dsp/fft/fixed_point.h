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
// policy, in the class docstring below and in docs/fft-design.md ("The
// fixed-point profiles").
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
//   J. W. Cooley, P. A. W. Lewis and P. D. Welch, "The fast Fourier
//     transform algorithm: Programming considerations in the calculation of
//     sine, cosine and Laplace transforms," J. Sound Vib. 12(3), 315-337,
//     1970 — the transform of 2M real samples through one M-point complex
//     transform: the real post-pass, the inverse's pre-pass and the
//     DC/Nyquist pair.
//   H. V. Sorensen, D. L. Jones, M. T. Heideman and C. S. Burrus,
//     "Real-valued fast Fourier transform algorithms," IEEE Trans. Acoust.,
//     Speech, Signal Process. ASSP-35(6), 849-863, 1987 — the same
//     half-length method and its inverse, among the real-valued algorithms.
//   The fixed-point arrangement of that method (which products, which
//   roundings, the coefficient table) is derived here, in real_post_pass's
//   docstring, from those equations and fft_arith.h; the complex kernel is
//   DspTap's own (Decision D2). How this pass was re-derived, under a
//   clean-room procedure (tap/DspTap#39): NOTICE.md and docs/fft-design.md,
//   "The fixed-point post-pass, re-derived".

#pragma once

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>
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
        /// a one-bit shift and the real post-pass (the half-length method of
        /// Cooley, Lewis and Welch 1970): per bin pair (k, M - k) one complex
        /// product by C_k = (1 + i W_N^k) / 2 (make_real_post_pass_table) and
        /// four exact sums, two roundings per output component; DC and
        /// Nyquist as the exact sum and difference of Z[0]'s two parts; bin
        /// N/4 untouched. The inverse runs the pre-pass first (the same
        /// statements with conj C_k, and DC/Nyquist halved with one one-bit
        /// rounding), then the conjugate kernel and the permutation.
        /// real_post_pass carries the derivation from the equations to the
        /// statements.
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
        /// (the golden model's unnormalized round trip has gain N/2).
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
        ///     returns e == 0. When e reaches the constant, the output is the
        ///     fixed schedule's bit for bit only if every stage shifted the
        ///     fixed amount (a block with no headroom at entry that stays at
        ///     full scale); a block that enters with headroom shifts less at
        ///     the first stage and catches up later, and then agrees with
        ///     the fixed output within the pinned bound (Q31 4 LSB at index
        ///     0, 2 LSB elsewhere; Q15 1 LSB, on a rounding tie of the
        ///     narrow), not bit for bit. `BfpAtTheFullExponentIsBitIdenticalToFixed`
        ///     (full-scale patterns), `BfpAtTheFullExponentAgreesWithFixedWithinPin`.
        ///
        /// Size: 4 <= N <= 65536, a power of two, fixed at construction
        /// (asserted, as for the floating profiles). Everything is allocated
        /// in the constructor: the three tables, and for Q15 an int32 work
        /// buffer of N (the in-place API is preserved at the caller's int16
        /// buffer; Q31 transforms in place with no work buffer). The
        /// transforms are noexcept and allocation-free; the object is
        /// copyable; there is no alignment requirement.
        ///
        /// Thread rule, as the engine trait k_is_shareable (Stage 4; a
        /// recorded deviation from audit Part 7, which lists both fixed
        /// profiles as shareable): Q31 true — a transform touches the
        /// caller's buffer and the const tables only, and its transforms are
        /// const so the compiler checks it; Q15 false — the widened block
        /// lives in the per-object work buffer m_work, and its transforms are
        /// non-const. Engine contract numbers otherwise: k_min_size = 4,
        /// k_max_size = 65536.
        ///
        /// Measured floors: the `[ floor ]` rows `NoiseFloorTracksWelchModel`
        /// prints (against basic_real_fft<double> on the same quantized
        /// input, output-referred, forward; N = 256 / 512 / 2048, white
        /// noise at 0 / -20 / -40 / -60 dBFS and an on-bin tone; the full
        /// table is in docs/fft-design.md beside the Welch model). Re-measured
        /// 2026-09-26 on the post-pass re-derived from the literature
        /// (real_post_pass): every number below is unchanged, the outputs
        /// being the same bit patterns (`OutputFingerprintIsPinned`), except
        /// the model ratio, whose post-pass term was re-derived with it:
        ///   - Q15, fixed: 0.26 - 0.30 LSB rms at every N and level, the
        ///     narrow's own round-half-up (1/sqrt(12) = 0.289); the int32
        ///     kernel's noise (0.58 LSB32 predicted) is 2^-14 of a Q15 LSB.
        ///     Per-bin SNR on full-scale white noise 69.1 / 66.5 / 60.2 dB at
        ///     N = 256 / 512 / 2048, falling 20 dB per 20 dB of input level:
        ///     the floor is level-independent, as Welch's model says. Largest
        ///     deviation over the adversarial sweep 0.50 LSB (pinned 1.0,
        ///     `SaturationFreeWorstCaseDoesNotWrap`).
        ///   - Q31, fixed: 0.67 - 0.75 LSB rms at every N and level. Welch's
        ///     variance-only model predicts 0.57 - 0.58 LSB; the measured /
        ///     model power ratio is 1.31 - 1.75 (pinned 3.51), the remainder
        ///     being the round-half-up bias of shr_round (+2^-(s+1) LSB per
        ///     s-bit shift: the one-bit shift before the post-pass puts +0.25
        ///     on every component entering it, which the pass carries to
        ///     between 0 and +0.5 per output component, bin by bin through
        ///     C_k; earlier stages' biases arrive with rotated phases).
        ///     Per-bin SNR on full-scale white noise 152.0 / 149.4 / 142.8 dB
        ///     at N = 256 / 512 / 2048. Largest deviation over the
        ///     adversarial sweep 4.25 LSB (inverse of a full-scale constant,
        ///     N = 2048, index 0; pinned 8.5).
        ///   - Block floating point: the floor follows the block, not the
        ///     format. Q15: 84 - 91 dB per-bin SNR on full-scale noise
        ///     (e = 5 - 7 where fixed takes 8 - 11), 78 - 86 dB at -40 dBFS
        ///     (e = 0); largest deviation 0.75 LSB (pinned 1.5). Q31: 157 -
        ///     161 dB at full scale, 154 - 156 dB at -40 dBFS, 135 - 136 dB
        ///     at -60 dBFS (e = 0: the input's own 2^-31 quantum is the
        ///     floor); 0.99 - 6.63 LSB rms as e falls from 9 to 0 (the
        ///     output is held at a larger scale; the SNR is the comparable
        ///     number); largest deviation 15.99 LSB (inverse of full-scale
        ///     binary noise, N = 1024, index 0; pinned 32 over the battery's
        ///     N = 4 / 8 / 16 / 64 / 512 / 1024 / 2048; above them the same
        ///     sweep measures 31.0 / 32.0 / 59.0 / 76.0 / 80.0 LSB at N =
        ///     4096 / 8192 / 16384 / 32768 / 65536, the DC-path bias below
        ///     growing with the exponent gap, pinned per size at 2x on the
        ///     host, `SaturationFreeWorstCaseIsPinnedPerLargeSize`).
        ///   - Honest limit: the round-half-up bias accumulates coherently
        ///     along the unrotated (DC) path, and under block floating point,
        ///     where a stage typically shifts one bit against a magnitude
        ///     gain of two, it grows about twice per stage instead of being
        ///     quartered. The battery sees it as the Q31 maxima above, all
        ///     at index 0 or 1, and as the mean of F(x) + F(-x)
        ///     (`RoundingBiasOnNegatedInputIsBounded`: 0.81 LSB Q31 fixed,
        ///     1.05 LSB Q31 block floating, pinned 1.63 / 2.1; the sum's
        ///     maximum 6 / 62 LSB, pinned 12 / 124, at index 0 / 1). Above
        ///     the sweep sizes the host-only
        ///     `RoundingBiasOnNegatedInputIsPinnedPerLargeSize` pins it per
        ///     size for N = 4096 ... 65536 at 2x the measured values (Q31
        ///     block floating: sum maximum 121 / 126 / 137 / 141 / 286 LSB,
        ///     mean 0.93 - 1.00 LSB; Q31 fixed 7 - 9 / 0.45 - 0.67 LSB;
        ///     docs/fft-design.md §3). Q15
        ///     never sees it (below the narrow's quantum). Convergent rounding
        ///     would remove it; that is a different trait contract
        ///     (fft_arith.h), not this kernel's choice.
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
            static constexpr bool        k_widens   = !std::is_same_v<Sample, wide>; ///< Q15: I/O width != kernel width
            static constexpr std::size_t k_min_size = 4;
            static constexpr std::size_t k_max_size = 65536;
            static constexpr bool        k_is_shareable = !k_widens; ///< Q31 true, Q15 false (class docstring)

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
                , m_twiddles(make_twiddle_table(n / 2))
                , m_post(make_real_post_pass_table(n))
                , m_bitrev(make_bit_reversal_table(n / 2))
                , m_work(std::is_same_v<Sample, std::int16_t> ? n : 0, wide{0}) {
                assert(n >= k_min_size && n <= k_max_size && (n & (n - 1)) == 0);
            }

            [[nodiscard]] std::size_t size() const noexcept { return m_n; }

            /// In-place forward transform: N samples -> packed spectrum, scaled
            /// by 2^-e. @return e (see the class docstring); nodiscard because
            /// under block floating point it is the output's scale. Q15 (the
            /// widening profile) runs through the per-object work buffer and
            /// is non-const; Q31 runs in place over the caller's buffer and the
            /// const tables, and is const — the k_is_shareable trait, checked
            /// by the compiler.
            [[nodiscard]] int forward_inplace(Sample* data) noexcept
                requires k_widens
            {
                return forward_on(enter(data), data);
            }
            [[nodiscard]] int forward_inplace(Sample* data) const noexcept
                requires(!k_widens)
            {
                return forward_on(data, data);
            }

            /// In-place inverse transform: packed spectrum -> N samples, the
            /// unnormalized inverse scaled by 2^-e. @return e. Constness as
            /// forward_inplace.
            [[nodiscard]] int inverse_inplace(Sample* data) noexcept
                requires k_widens
            {
                return inverse_on(enter(data), data);
            }
            [[nodiscard]] int inverse_inplace(Sample* data) const noexcept
                requires(!k_widens)
            {
                return inverse_on(data, data);
            }

          private:
            /// The forward transform over the kernel's view `a` of the caller's
            /// buffer `data` (the work buffer for Q15, `data` itself for Q31).
            int forward_on(wide* a, Sample* data) const noexcept {
                // The fixed policy's cumulative shift is tracked beside the
                // exponent actually applied; the two coincide under
                // scaling::fixed and the first bounds the second under
                // scaling::block_floating (stage_shift).
                int e   = 0;
                int cum = arith::k_fixed_scaling_input_pre_shift; // folded into the first stage's shift
                e       = complex_kernel<false>(a, e, cum);
                permute(a);
                cum += 1; // the real post-pass: one bit
                e += shift_block(a, stage_shift(a, 1, e, cum));
                real_post_pass<false>(a);
                return leave(a, data, e, cum);
            }

            /// The inverse transform over the kernel's view, as forward_on.
            int inverse_on(wide* a, Sample* data) const noexcept {
                int e   = 0;
                int cum = arith::k_fixed_scaling_input_pre_shift + 1; // pre-shift and the pre-pass bit
                e += shift_block(a, stage_shift(a, 1, e, cum));
                real_post_pass<true>(a);
                e = complex_kernel<true>(a, e, cum);
                permute(a);
                return leave(a, data, e, cum);
            }

            /// Q15 only: widen the caller's block into the work buffer, the
            /// kernel's view of it (Q31's view is the caller's buffer itself,
            /// passed straight through by the const overloads above).
            wide* enter(Sample* data) noexcept
                requires k_widens
            {
                for (std::size_t i = 0; i < m_n; ++i) {
                    m_work[i] = arith::widen(data[i]);
                }
                return m_work.data();
            }

            /// Hands the block back to the caller. Q15 narrows Q2.29 -> Q0.15
            /// (round-half-up, saturating); under block floating point the
            /// block is first shifted so that it fits below full scale, one
            /// more stage of growth 1 whose shift is 0 under scaling::fixed
            /// (the fixed schedule already put the result within +-1.0).
            int leave(wide* a, Sample* data, int e, int cum) const noexcept {
                if constexpr (k_widens) {
                    e += shift_block(a, stage_shift(a, 1, e, cum));
                    for (std::size_t i = 0; i < m_n; ++i) {
                        data[i] = arith::narrow(a[i]);
                    }
                }
                else {
                    (void)a;
                    (void)data;
                    (void)cum;
                }
                return e;
            }

            /// The shift applied before a stage whose outputs can grow by
            /// `growth` bits, given the exponent applied so far and the fixed
            /// policy's cumulative shift through this stage.
            ///   scaling::fixed:          cum - e, i.e. exactly the stage's growth
            ///                            (plus the input pre-shift on the first).
            ///   scaling::block_floating: the headroom scan's demand,
            ///                            growth + 1 - headroom (one bit is never
            ///                            consumed: fft_arith.h), clamped to
            ///                            [0, cum - e]. The upper clamp is what
            ///                            makes e <= fixed_scaling_exponent(N)
            ///                            a contract; it binds only when the
            ///                            block has caught up with the fixed
            ///                            schedule (e == cum before this stage)
            ///                            and reads h = 0, where the fixed
            ///                            schedule's saturation-freedom proof
            ///                            applies. Note that a block can be
            ///                            BEHIND the fixed schedule (an input
            ///                            with headroom at entry: h = 1 shifts
            ///                            2 where fixed shifts 3) and catch up
            ///                            with a larger shift later, so
            ///                            reaching the constant does not mean
            ///                            the fixed rounding history was
            ///                            followed (class docstring).
            int stage_shift(const wide* a, int growth, int e, int cum) const noexcept {
                const int cap = cum - e;
                if constexpr (k_block_floating) {
                    const int headroom = work::headroom_bits(a, m_n);
                    return std::clamp(growth + 1 - headroom, 0, cap);
                }
                else {
                    (void)a;
                    (void)growth;
                    return cap;
                }
            }

            /// shr_round every value of the block by `bits`; returns bits.
            int shift_block(wide* a, int bits) const noexcept {
                if (bits > 0) {
                    for (std::size_t i = 0; i < m_n; ++i) {
                        a[i] = work::shr_round(a[i], bits);
                    }
                }
                return bits;
            }

            /// (xr + i xi) * (wr + i wi), or by the conjugate twiddle when
            /// Inverse: two mul_coeff roundings per output component, exactly
            /// the form fft_arith.h specifies (no fused variant).
            template <bool Inverse>
            static void rotate(wide& xr, wide& xi, coeff wr, coeff wi) noexcept {
                const wide rr = work::mul_coeff(xr, wr);
                const wide ii = work::mul_coeff(xi, wi);
                const wide ri = work::mul_coeff(xr, wi);
                const wide ir = work::mul_coeff(xi, wr);
                if constexpr (Inverse) {
                    xr = work::add(rr, ii);
                    xi = work::sub(ir, ri);
                }
                else {
                    xr = work::sub(rr, ii);
                    xi = work::add(ri, ir);
                }
            }

            /// The radix-4 decimation-in-frequency kernel of length M = N/2 on
            /// the interleaved block, radix-2 final stage when log2 M is odd,
            /// output in bit-reversed order (permute() follows). Each stage
            /// shifts its inputs by stage_shift() before the butterfly (2 bits
            /// of growth per radix-4 stage, 1 per radix-2). Returns e.
            template <bool Inverse>
            int complex_kernel(wide* a, int e, int& cum) const noexcept {
                const std::size_t m = m_n / 2;
                std::size_t       span;
                for (span = m; span >= 4; span /= 4) {
                    cum += 2;
                    const int s = stage_shift(a, 2, e, cum);
                    e += s;
                    radix4_stage<Inverse>(a, span, m / span, s);
                }
                if (span == 2) {
                    cum += 1;
                    const int s = stage_shift(a, 1, e, cum);
                    e += s;
                    radix2_stage(a, s);
                }
                return e;
            }

            /// One radix-4 DIF stage over sub-transforms of `span` complex
            /// values. For each butterfly (inputs at q, q + L/4, q + L/2,
            /// q + 3L/4, each shr_round'ed by `shift` first), with
            /// t0 = a0 + a2, t1 = a0 - a2, t2 = a1 + a3, t3 = a1 - a3:
            ///   y0 = t0 + t2               -> q
            ///   y1 = (t1 + i t3) W^q        -> q + L/2   (W_4 = +i forward, -i inverse)
            ///   y2 = (t0 - t2) W^2q         -> q + L/4
            ///   y3 = (t1 - i t3) W^3q       -> q + 3L/4
            /// y1 and y2 swap places relative to the textbook radix-4 so the
            /// final order is a plain bit reversal (radix-2^2 placement),
            /// W = exp(+2*pi*i/L) forward and its conjugate inverse. q = 0
            /// skips the rotation: mul_coeff by k_coeff_one is the identity,
            /// so the skip is bit-identical. `stride` is M/L, the table step.
            template <bool Inverse>
            void radix4_stage(wide* a, std::size_t span, std::size_t stride, int shift) const noexcept {
                const std::size_t quarter = span / 4;
                const coeff*      w       = m_twiddles.data();
                for (std::size_t block = 0; block < m_n / 2; block += span) {
                    for (std::size_t q = 0; q < quarter; ++q) {
                        wide* const p0 = a + 2 * (block + q);
                        wide* const p1 = p0 + 2 * quarter;
                        wide* const p2 = p1 + 2 * quarter;
                        wide* const p3 = p2 + 2 * quarter;

                        const wide a0r = work::shr_round(p0[0], shift);
                        const wide a0i = work::shr_round(p0[1], shift);
                        const wide a1r = work::shr_round(p1[0], shift);
                        const wide a1i = work::shr_round(p1[1], shift);
                        const wide a2r = work::shr_round(p2[0], shift);
                        const wide a2i = work::shr_round(p2[1], shift);
                        const wide a3r = work::shr_round(p3[0], shift);
                        const wide a3i = work::shr_round(p3[1], shift);

                        const wide t0r = work::add(a0r, a2r);
                        const wide t0i = work::add(a0i, a2i);
                        const wide t1r = work::sub(a0r, a2r);
                        const wide t1i = work::sub(a0i, a2i);
                        const wide t2r = work::add(a1r, a3r);
                        const wide t2i = work::add(a1i, a3i);
                        const wide t3r = work::sub(a1r, a3r);
                        const wide t3i = work::sub(a1i, a3i);

                        // y1 = t1 + i t3 (forward) / t1 - i t3 (inverse); y3 the other.
                        wide y1r;
                        wide y1i;
                        wide y3r;
                        wide y3i;
                        if constexpr (Inverse) {
                            y1r = work::add(t1r, t3i);
                            y1i = work::sub(t1i, t3r);
                            y3r = work::sub(t1r, t3i);
                            y3i = work::add(t1i, t3r);
                        }
                        else {
                            y1r = work::sub(t1r, t3i);
                            y1i = work::add(t1i, t3r);
                            y3r = work::add(t1r, t3i);
                            y3i = work::sub(t1i, t3r);
                        }
                        wide y2r = work::sub(t0r, t2r);
                        wide y2i = work::sub(t0i, t2i);

                        if (q != 0) {
                            const std::size_t k1 = 2 * (q * stride);
                            const std::size_t k2 = 2 * k1;
                            const std::size_t k3 = 3 * k1;
                            rotate<Inverse>(y1r, y1i, w[k1], w[k1 + 1]);
                            rotate<Inverse>(y2r, y2i, w[k2], w[k2 + 1]);
                            rotate<Inverse>(y3r, y3i, w[k3], w[k3 + 1]);
                        }

                        p0[0] = work::add(t0r, t2r);
                        p0[1] = work::add(t0i, t2i);
                        p2[0] = y1r; // X[4m+1] sequence at q + L/2
                        p2[1] = y1i;
                        p1[0] = y2r; // X[4m+2] sequence at q + L/4
                        p1[1] = y2i;
                        p3[0] = y3r;
                        p3[1] = y3i;
                    }
                }
            }

            /// The radix-2 final stage (span 2, twiddle exactly 1): y0 = a0 + a1
            /// at q, y1 = a0 - a1 at q + 1, inputs shr_round'ed by `shift`.
            /// Direction-free: with W_2 = -1 there is nothing to conjugate.
            void radix2_stage(wide* a, int shift) const noexcept {
                for (std::size_t i = 0; i < m_n; i += 4) {
                    const wide a0r = work::shr_round(a[i], shift);
                    const wide a0i = work::shr_round(a[i + 1], shift);
                    const wide a1r = work::shr_round(a[i + 2], shift);
                    const wide a1i = work::shr_round(a[i + 3], shift);
                    a[i]           = work::add(a0r, a1r);
                    a[i + 1]       = work::add(a0i, a1i);
                    a[i + 2]       = work::sub(a0r, a1r);
                    a[i + 3]       = work::sub(a0i, a1i);
                }
            }

            /// Bit-reversal permutation of the M complex values (no arithmetic).
            void permute(wide* a) const noexcept {
                const std::size_t m = m_n / 2;
                for (std::size_t i = 0; i < m; ++i) {
                    const std::size_t j = m_bitrev[i];
                    if (i < j) {
                        std::swap(a[2 * i], a[2 * j]);
                        std::swap(a[2 * i + 1], a[2 * j + 1]);
                    }
                }
            }

            /// The real post-pass (forward: after the kernel, the permutation and
            /// the one-bit shift) and pre-pass (Inverse: after the one-bit shift,
            /// before the conjugate kernel). The half-length method of Cooley,
            /// Lewis and Welch (1970) and Sorensen et al. (1987), in the
            /// library's convention and arranged for fft_arith's operations.
            ///
            /// Forward. With M = N/2, z[n] = x[2n] + i x[2n+1] and Z the kernel's
            /// length-M transform (W_M = exp(+2*pi*i/M)), the even and odd
            /// samples' transforms are E[k] = (Z[k] + conj Z[M-k]) / 2 and
            /// O[k] = (Z[k] - conj Z[M-k]) / (2i), both conjugate-symmetric,
            /// and X[k] = E[k] + W_N^k O[k], X[M-k] = conj(E[k] - W_N^k O[k]).
            /// The block entering here is a[k] = Z[k] * 2^-(s+1) (s the kernel's
            /// shift, 1 the pass's own bit, which pays for its growth); the pass
            /// must leave Y[k] = X[k] * 2^-(s+1), i.e.
            /// Y[k] = (u + v) / 2 + (-i W_N^k)(u - v) / 2 with u = a[k],
            /// v = conj a[M-k]. Writing (u + v) / 2 as u - (u - v) / 2 folds
            /// the 1/2 and the rotation into ONE complex product:
            ///   G = C_k (u - v),  C_k = (1 + i W_N^k) / 2   (make_real_post_pass_table)
            ///   Y[k]   = u - G                   (= E + W O)
            ///   Y[M-k] = a[M-k] + conj G         (= conj(v + G) = conj(E - W O))
            /// per pair (k, M-k), 1 <= k < M/2: one sub and one add to form
            /// u - v, the two-rounding complex multiply (rotate), and four exact
            /// add/sub. Each output component carries TWO roundings (the two
            /// mul_coeff of G) on top of the one-bit shift before the pass; no
            /// shr_round is needed for the 1/2 factors, which live in C_k.
            ///   - DC and Nyquist (k = 0): Z[0]'s partner is itself, E[0] =
            ///     Re Z[0], O[0] = Im Z[0], W_N^0 = 1, W_N^M = -1, so
            ///     X[0] = Re Z[0] + Im Z[0] and X[M] = Re Z[0] - Im Z[0]; the block
            ///     already carries the 2^-(s+1), so Y[0] = a[0] + a[1] and
            ///     Y[M] = a[0] - a[1] (packed at a[0], a[1]): exact, no rounding.
            ///   - Bin N/4 (k = M/2, its own partner): W_N^(N/4) = i, C = 0, so
            ///     Y[M/2] = a[M/2]: untouched, exact.
            ///
            /// Inverse. Given X, Z[k] = E[k] + i O[k] with E[k] = (X[k] +
            /// conj X[M-k]) / 2 and W_N^k O[k] = (X[k] - conj X[M-k]) / 2 makes
            /// the unnormalized length-M inverse kernel return
            /// sum_k Z[k] W_M^-nk = M z[n] = (N/2)(x[2n] + i x[2n+1]), which is
            /// the golden model's unnormalized inverse; so the pre-pass
            /// carries no extra factor: entering a[k] = X[k] * 2^-(p+1) (p the
            /// input pre-shift, 1 the pass's bit), it leaves Z[k] * 2^-(p+1).
            /// With u = a[k], v = conj a[M-k]:
            /// Z[k] = (u + v) / 2 + i W_N^-k (u - v) / 2 = u - conj(C_k)(u - v),
            /// and conj Z[M-k] = v + conj(C_k)(u - v): the forward's statements
            /// with the conjugate coefficient (rotate<true>), the same table.
            ///   - k = 0: Z[0] = ((X[0] + X[M]) / 2, (X[0] - X[M]) / 2), so
            ///     a[0] <- shr_round(a[0] + a[1], 1), a[1] <- shr_round(a[0] - a[1], 1):
            ///     ONE one-bit rounding per component (the inverse's DC/Nyquist
            ///     has a 1/2 that no coefficient absorbs).
            ///   - k = M/2: conj(C) = 0, a[M/2] untouched.
            ///
            /// Bounds (docs/fft-design.md §3 is the proof): the entering block
            /// has complex magnitude <= 2^29.5 (plus rounding slack) under both
            /// policies; u - v <= 2^30.5, the partial products <= 2^29.5,
            /// |G| <= 2^30 (|C_k| <= 1/sqrt(2)), and every output is <= 2^30 in
            /// magnitude because the pair is energy-preserving,
            /// |Y[k]|^2 + |Y[M-k]|^2 = (|u + v|^2 + |u - v|^2) / 2 = |u|^2 + |v|^2
            /// (the same for the pre-pass); DC / Nyquist sums <= 2^30, the
            /// inverse's halved ones <= 2^29.
            template <bool Inverse>
            void real_post_pass(wide* a) const noexcept {
                const std::size_t m = m_n / 2;
                const coeff*      c = m_post.data();
                const wide        r = a[0];
                const wide        i = a[1];
                if constexpr (Inverse) {
                    a[0] = work::shr_round(work::add(r, i), 1);
                    a[1] = work::shr_round(work::sub(r, i), 1);
                }
                else {
                    a[0] = work::add(r, i);
                    a[1] = work::sub(r, i);
                }
                for (std::size_t k = 1; k < m / 2; ++k) {
                    wide* const pk = a + 2 * k;               // u = a[k]
                    wide* const pj = a + 2 * (m - k);         // v = conj a[M - k]
                    wide        gr = work::sub(pk[0], pj[0]); // (u - v), real
                    wide        gi = work::add(pk[1], pj[1]); // (u - v), imaginary
                    rotate<Inverse>(gr, gi, c[2 * k], c[2 * k + 1]);
                    pk[0] = work::sub(pk[0], gr);
                    pk[1] = work::sub(pk[1], gi);
                    pj[0] = work::add(pj[0], gr);
                    pj[1] = work::sub(pj[1], gi);
                }
            }

            static constexpr int log2_of(std::size_t n) noexcept { return std::bit_width(n) - 1; }

            std::size_t                m_n;
            std::vector<coeff>         m_twiddles; ///< W_M^k, M = N/2, interleaved (tables.h)
            std::vector<coeff>         m_post;     ///< C_k = (1 + i W_N^k) / 2, k < N/4, interleaved (tables.h)
            std::vector<std::uint32_t> m_bitrev;   ///< permutation of the M complex outputs
            std::vector<wide>          m_work;     ///< Q15 only: the widened block
        };

    } // namespace detail

} // namespace tap::dsp
