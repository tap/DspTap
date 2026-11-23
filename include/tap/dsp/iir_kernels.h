/// @file iir_kernels.h
/// @brief Recursive kernels per sample format: allpass sections, the biquad cascade, the two-path and four-path steps.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The kernel half of SampleRateTap's recursive polyphase stages
// (docs/IIR_PLAN.md, I1, I8, I9, I12, I14), over the designs of
// iir_design.h, in the four formats of sample_traits.h: double (the golden
// model), float (the embedded floating profile), Q15 and Q31. Published
// literature only (iir_design.h's sources; Jackson 1969 and 1970 for the
// fixed-point rounding and scaling).
//
// Where the line is drawn (fir_kernels.h's line): a kernel computes outputs
// from the inputs its paths take for them, against state the caller owns
// (dot_row's history window is the FIR analogue), and keeps nothing of its
// own. The commutator schedule — which input goes to which path and when, the
// one-sample odd latch of a decimator (x[2k - 1]), the four-path's three, the
// near-linear delay path's history (x[2(k - K)]) — the stream's phase across
// calls, channel interleaving, outputs_for, flush and the latency members are
// the stage machine's (rational's basic_iir_stage, IIR-b), exactly as a FIR
// stage keeps its history and phase and calls dot_row. The recursive states
// themselves are part of the structure's arithmetic (their headroom and
// rounding are the format contract below), so the kernels update them; the
// caller only stores them, per channel, and zeroes them on reset.
//
// The structures (iir_design.h): one-multiplier first-order allpass sections
// (a + w^-1) / (1 + a w^-1), y = x1 + a (x - y1), where a cascade shares
// states (section i's previous output is section i+1's previous input, so n
// sections hold n + 1 states); second-order allpass sections (c2 + c1 w^-1 +
// w^-2) / (1 + c1 w^-1 + c2 w^-2), y = c2 (x - y2) + c1 (x1 - y1) + x2, two
// states per joint, shared the same way; biquads (b0 + b1 z^-1 + b2 z^-2) /
// (1 + a1 z^-1 + a2 z^-2), a first-order section being one with b2 = a2 = 0.
//
// Contract, as numbers (iir_traits<S> carries them as constants):
//   - double: states, coefficients and products in double; the golden model.
//   - float: states, coefficients and products in float (I14, a documented
//     exception to DspTap's "fragile recursions run in double" rule: the
//     one-multiplier section is structurally passive and insensitive to its
//     coefficient, and IIR-0 measured the worst error with float states at
//     -137.1 dB re full scale, the near-linear 120 dB, 5/12 interpolator on
//     white noise, so 120 dB stays offered in float; float states count 28
//     times fewer instructions than double on the soft-double M33 and M4F).
//     Biquads in transposed direct form II.
//   - Q15 (int16_t samples, Q0.15): states int32 in Q2.29 (two integer guard
//     bits, |state| < 4: the stored states' L1 norms reach 2.69 and 3.71 at
//     the 70 dB and 120 dB spec points, IIR-0); the input enters as x << 14.
//     One-multiplier coefficients int32 in Q0.31 (every a in [0, 1));
//     second-order allpass and biquad coefficients int32 in Q2.30 (|c| < 2).
//     Every product and sum in 64 bits; the difference nodes x - y1, x - y2,
//     x1 - y1 are formed there and never stored (they reach 5.22 in L1 at
//     120 dB, past two guard bits: IIR-0). Biquads in direct form I with
//     shared joints. One rounding point per section, at its output: floor
//     truncation (an arithmetic shift) in the one-multiplier sections, which
//     leaves no zero-input limit cycle for a >= 0 (its recursion u <-
//     floor(-a u) alternates sign and shrinks; IIR-0 searched every section
//     exhaustively); magnitude truncation (toward zero) in every
//     second-order section, allpass or biquad, where floor truncation leaves
//     limit cycles in 28 of 88 and magnitude truncation in none (IIR-0, by
//     Jackson's 1969 bound). Round to nearest is never used inside the
//     recursion: it cycles in 44 of 74 first-order and 65 of 88 second-order
//     sections. The output leaves Q2.29 once, round half up, saturating to
//     Q15: the path output >> 14, the two-path half-sum >> 15, the four-path
//     quarter-sum >> 16.
//   - Q31 (int32_t samples, Q0.31): states int64 in Q2.60 (two guard bits;
//     the 64-bit word's top bit is the difference node's), the input as
//     x << 29; coefficients as for Q15 (Q0.31 and Q2.30 in int32). The
//     products are 64 x 32 = 96-bit, formed exactly from two 32 x 32
//     multiplies, summed exactly, rounded once (the same truncations as Q15),
//     so the result is the exact product rounded at the section's output.
//     The output leaves Q2.60 once: >> 29 (and >> 30, >> 31 for the half- and
//     quarter-sums), round half up, saturating to Q31.
//   - Saturation happens at the output only. Inside, the headroom is the
//     design's: a full-scale input (square wave, clipped sine) never takes a
//     state past |4| at the tabulated designs (tests/test_iir_fixed.cpp pins
//     it); an input that did would wrap, modulo 2^32 or 2^64, never invoke
//     undefined behaviour (sums run in unsigned arithmetic).
//   - The floors (tests/test_iir_fixed.cpp): the Q15 output equals the double
//     structure rounded to Q15 within one LSB, the recursion adding nothing
//     measurable above the output's own rounding (IIR-0: -98.1 dB for a
//     -1 dBFS sine at 70 dB, 3/8).
//   - Coefficient tables in each format (nonlinear_halfband_table,
//     near_linear_halfband_table, fourpath_table, biquad_cascade_table) are
//     built from the double designs by constexpr constructors, so a named
//     design's fixed-point tables round from the same compile-time doubles
//     (I4); they own fixed-capacity arrays (no allocation) and state the
//     per-channel state size their kernels read and write.
//
// Real-time safe: every kernel is noexcept, allocation-free and branch-light
// over caller-owned memory. The structure kernels (halfband_decimate,
// halfband_interpolate, near_linear_*, fourpath_*, biquad_cascade) state the
// flush-to-zero precondition (subnormals.h: TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED,
// a debug-build read of the control register where a bit governs S, never a
// write); the section and path primitives do not check.
//
// Honest limits (pinned by tests/test_iir_kernels.cpp):
//   - Non-finite input poisons a recursive state until the caller zeroes it
//     (an FIR stage clears a non-finite burst after one window; a recursion
//     never does). I9 decided against a per-sample check.
//   - Subnormals: without the host's flush-to-zero, or on a target with no
//     bit for the type (subnormals.h's table), a state decays to a residue of
//     a few multiples of the smallest subnormal and stays there (measured: 2
//     at 70 dB 3/8, 5 to 7 at 120 dB 5/12 for the half-bands; 4 to 20 for the
//     Chebyshev II cascades), and on x86 the silence after signal runs about
//     ten times slower. Under flush-to-zero the allpass sections (both
//     half-band variants, the four-path stage) reach an exactly zero state,
//     with one exception: the floating near-linear half-band at 120 dB, 5/12,
//     divisor 1 (N19, nine second-order sections) holds a zero-input limit
//     cycle under flush-to-zero, at 3.7 to 4.1 times the smallest normal in
//     float and double, from 12 to 20 of 39 noise starts (SampleRateTap's
//     IIR-b; rational's iir_stage.h states it and pins its bound): a
//     second-order allpass section's product flushes and the section lifts
//     the state again, the mechanism of the next item.
//   - The floating biquad cascade under flush-to-zero: once |a2 y| falls
//     below the smallest normal its product flushes, and for that sample the
//     section runs as a first-order recursion with |a1| > 1, which lifts it
//     again; so the guard band filter's Chebyshev II of order 7 or more can
//     hold a zero-input limit cycle at a few times the smallest normal and
//     never reach an exactly zero state, from some states and not others:
//     measured from full-scale square-wave starts, up to 13.1 times in this
//     transposed direct form II (about 3e-37 in float, 3e-307 in double),
//     and direct form I in floating point does the same (up to 7.1 times).
//     Orders 3 and 5 at 1/3 drain to zero. The fixed-point cascade
//     (magnitude truncation) always drains. A guard band stage's flush
//     window (IIR_PLAN I6) therefore cannot promise an exactly zero state
//     in floating point at those orders. Decided (SampleRateTap's
//     docs/IIR_PLAN.md question 19): an honest limit, kept, with no deadband
//     in these kernels and no design dropped. A consumer states the residue
//     its stage leaves after its flush window and pins a bound on it, as
//     rational's guard band filter does at order 7 (its guard_band.h);
//     DspTap's own measurement is the one above
//     (tests/test_iir_kernels.cpp).
//   - Contraction: the floating kernels are written as plain expressions, so
//     a compiler that contracts them into FMAs (GCC by default where the
//     target has one, Clang within a statement) moves float and double
//     output bits per target; their contracts are tolerances against the
//     double golden model, never bits. The fixed-point outputs are integers
//     and bit-identical everywhere (tests/test_iir_fixed.cpp pins them).

#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>

#include "tap/dsp/detail/expects.h"
#include "tap/dsp/iir_design.h"
#include "tap/dsp/sample_traits.h"
#include "tap/dsp/subnormals.h"

// The no-alias qualifier fir_kernels.h defines, spelled identically so the
// two headers may be included in either order.
#if !defined(TAP_DSP_RESTRICT)
#if defined(_MSC_VER)
#define TAP_DSP_RESTRICT __restrict
#else
#define TAP_DSP_RESTRICT __restrict__
#endif
#endif

/// Keeps biquad_cascade_block's group loop out of line (its comment says why).
#if !defined(TAP_DSP_IIR_NOINLINE)
#if defined(__GNUC__) || defined(__clang__)
#define TAP_DSP_IIR_NOINLINE [[gnu::noinline]]
#elif defined(_MSC_VER)
#define TAP_DSP_IIR_NOINLINE __declspec(noinline)
#else
#define TAP_DSP_IIR_NOINLINE
#endif
#endif

namespace tap::dsp {

    namespace detail {

        /// v / 2^shift rounded toward zero (magnitude truncation).
        constexpr std::int64_t shift_toward_zero(std::int64_t v, int shift) noexcept {
            const std::int64_t bias = (v >> 63) & ((std::int64_t{1} << shift) - 1);
            return (v + bias) >> shift;
        }

        /// a + b modulo 2^64: exact whenever the true sum fits int64, and
        /// never undefined.
        constexpr std::int64_t wrap_add(std::int64_t a, std::int64_t b) noexcept {
            return static_cast<std::int64_t>(static_cast<std::uint64_t>(a) + static_cast<std::uint64_t>(b));
        }

        /// An exact 96-bit accumulator of int32 x int64 products: the value is
        /// hi 2^32 + lo, 0 <= lo < 2^35 until normalized.
        struct wide_accumulator {
            std::uint64_t hi = 0; ///< two's complement, modulo 2^64
            std::uint64_t lo = 0;

            /// += c d, exactly: d = dh 2^32 + dl with dl in [0, 2^32).
            constexpr void mac(std::int32_t c, std::int64_t d) noexcept {
                const std::int64_t dh = d >> 32;
                const std::int64_t dl = d & 0xffffffff;
                const std::int64_t ph = static_cast<std::int64_t>(c) * dh; // |.| < 2^62
                const std::int64_t pl = static_cast<std::int64_t>(c) * dl; // |.| < 2^63
                hi += static_cast<std::uint64_t>(ph) + static_cast<std::uint64_t>(pl >> 32);
                lo += static_cast<std::uint64_t>(pl & 0xffffffff);
            }

            /// += v 2^30, exactly.
            constexpr void add_shifted30(std::int64_t v) noexcept {
                hi += static_cast<std::uint64_t>(v >> 2);
                lo += static_cast<std::uint64_t>(v & 3) << 30;
            }

            /// The value / 2^shift for shift in [1, 32]: floor, or toward zero.
            template <bool TowardZero>
            constexpr std::int64_t shifted(int shift) const noexcept {
                const std::uint64_t h    = hi + (lo >> 32);
                const std::uint64_t l    = lo & 0xffffffff;
                const auto          sh   = static_cast<std::int64_t>(h);
                const std::uint64_t mask = (std::uint64_t{1} << shift) - 1;
                const auto          fl   = static_cast<std::int64_t>((h << (32 - shift)) + (l >> shift));
                if constexpr (TowardZero) {
                    return fl + ((sh < 0 && (l & mask) != 0) ? 1 : 0);
                }
                else {
                    return fl;
                }
            }
        };

        /// floor(a d / 2^31) for an int32 a and an int64 d, exactly.
        constexpr std::int64_t mul_q31_floor(std::int32_t a, std::int64_t d) noexcept {
            const std::int64_t dh = d >> 32;
            const std::int64_t dl = d & 0xffffffff;
            const std::int64_t ph = static_cast<std::int64_t>(a) * dh; // |.| < 2^62
            const std::int64_t pl = static_cast<std::int64_t>(a) * dl; // |.| < 2^63
            // a d = ph 2^32 + pl, so a d / 2^31 = 2 ph + pl / 2^31.
            return wrap_add(static_cast<std::int64_t>(static_cast<std::uint64_t>(ph) << 1), pl >> 31);
        }

    } // namespace detail

    /// Primary template intentionally undefined: the four formats below.
    template <sample_type S>
    struct iir_traits;

    namespace detail {

        /// The floating formats: everything in F (I14 for float).
        template <typename F>
        struct floating_iir_traits {
            using state  = F; ///< path and section values
            using coeff1 = F; ///< one-multiplier allpass coefficient
            using coeff2 = F; ///< second-order allpass and biquad coefficient

            static constexpr bool k_is_fixed_point = false;

            static constexpr coeff1 make_coeff1(double a) noexcept { return static_cast<F>(a); }
            static constexpr coeff2 make_coeff2(double c) noexcept { return static_cast<F>(c); }
            static constexpr state  input(F x) noexcept { return x; }
            static constexpr F      output(state y) noexcept { return y; }
            static constexpr F      output_mean2(state a, state b) noexcept { return F(0.5) * (a + b); }
            static constexpr F      output_mean4(state a, state b, state c, state d) noexcept {
                return F(0.25) * ((a + b) + (c + d));
            }

            /// y = x1 + a (x - y1).
            static constexpr state allpass1(state x, state x1, state y1, coeff1 a) noexcept {
                return x1 + a * (x - y1);
            }

            /// y = c2 (x - y2) + c1 (x1 - y1) + x2.
            static constexpr state allpass2(state x, state x1, state x2, state y1, state y2, coeff2 c1,
                                            coeff2 c2) noexcept {
                return c2 * (x - y2) + c1 * (x1 - y1) + x2;
            }
        };

    } // namespace detail

    /// double: the golden model.
    template <>
    struct iir_traits<double> : detail::floating_iir_traits<double> {};

    /// float: float states and products (I14).
    template <>
    struct iir_traits<float> : detail::floating_iir_traits<float> {};

    // ANCHOR: iir_q15
    /// Q15: Q0.15 samples, Q2.29 int32 states, Q0.31 / Q2.30 coefficients,
    /// 64-bit products and sums, one rounding per section (the file header).
    template <>
    struct iir_traits<std::int16_t> {
        using state  = std::int32_t;
        using coeff1 = std::int32_t;
        using coeff2 = std::int32_t;

        static constexpr bool k_is_fixed_point   = true;
        static constexpr int  k_state_frac_bits  = 29; ///< Q2.29: two integer guard bits
        static constexpr int  k_input_shift      = 14; ///< Q0.15 -> Q2.29
        static constexpr int  k_coeff1_frac_bits = 31; ///< Q0.31 one-multiplier coefficients
        static constexpr int  k_coeff2_frac_bits = 30; ///< Q2.30 second-order coefficients
        static_assert(k_state_frac_bits - k_input_shift == 15, "Q15 ladder: Q0.15 << 14 = Q2.29");

        static constexpr coeff1 make_coeff1(double a) noexcept {
            return detail::round_sat<std::int32_t>(a * static_cast<double>(std::int64_t{1} << k_coeff1_frac_bits));
        }
        static constexpr coeff2 make_coeff2(double c) noexcept {
            return detail::round_sat<std::int32_t>(c * static_cast<double>(std::int64_t{1} << k_coeff2_frac_bits));
        }
        static constexpr state input(std::int16_t x) noexcept {
            return static_cast<std::int32_t>(x) * (std::int32_t{1} << k_input_shift);
        }
        /// Q2.29 -> Q15, round half up, saturating.
        static constexpr std::int16_t output(state y) noexcept {
            return detail::clamp_sat<std::int16_t>((static_cast<std::int64_t>(y) + (std::int64_t{1} << 13)) >> 14);
        }
        static constexpr std::int16_t output_mean2(state a, state b) noexcept {
            const std::int64_t s = static_cast<std::int64_t>(a) + b;
            return detail::clamp_sat<std::int16_t>((s + (std::int64_t{1} << 14)) >> 15);
        }
        static constexpr std::int16_t output_mean4(state a, state b, state c, state d) noexcept {
            const std::int64_t s = (static_cast<std::int64_t>(a) + b) + (static_cast<std::int64_t>(c) + d);
            return detail::clamp_sat<std::int16_t>((s + (std::int64_t{1} << 15)) >> 16);
        }

        /// y = x1 + floor(a (x - y1) / 2^31): the difference in 64 bits, the
        /// Q0.31 x Q3.29 product below 2^63, floor truncation.
        static constexpr state allpass1(state x, state x1, state y1, coeff1 a) noexcept {
            const std::int64_t p = static_cast<std::int64_t>(a) * (static_cast<std::int64_t>(x) - y1);
            return static_cast<std::int32_t>(static_cast<std::int64_t>(x1) + (p >> 31));
        }

        /// y = trunc((x2 2^30 + c2 (x - y2) + c1 (x1 - y1)) / 2^30), toward zero.
        static constexpr state allpass2(state x, state x1, state x2, state y1, state y2, coeff2 c1,
                                        coeff2 c2) noexcept {
            std::int64_t full = static_cast<std::int64_t>(x2) * (std::int64_t{1} << 30);
            full = detail::wrap_add(full, static_cast<std::int64_t>(c2) * (static_cast<std::int64_t>(x) - y2));
            full = detail::wrap_add(full, static_cast<std::int64_t>(c1) * (static_cast<std::int64_t>(x1) - y1));
            return static_cast<std::int32_t>(detail::shift_toward_zero(full, 30));
        }

        /// Direct form I: y = trunc((b0 x + b1 x1 + b2 x2 + na1 y1 + na2 y2) / 2^30),
        /// na = -a, toward zero.
        static constexpr state biquad(state x, state x1, state x2, state y1, state y2,
                                      const std::array<coeff2, 5>& c) noexcept {
            std::int64_t full = static_cast<std::int64_t>(c[0]) * x;
            full              = detail::wrap_add(full, static_cast<std::int64_t>(c[1]) * x1);
            full              = detail::wrap_add(full, static_cast<std::int64_t>(c[2]) * x2);
            full              = detail::wrap_add(full, static_cast<std::int64_t>(c[3]) * y1);
            full              = detail::wrap_add(full, static_cast<std::int64_t>(c[4]) * y2);
            return static_cast<std::int32_t>(detail::shift_toward_zero(full, 30));
        }
    };
    // ANCHOR_END: iir_q15

    // ANCHOR: iir_q31
    /// Q31: Q0.31 samples, Q2.60 int64 states, Q0.31 / Q2.30 coefficients,
    /// exact 96-bit products, one rounding per section (the file header).
    template <>
    struct iir_traits<std::int32_t> {
        using state  = std::int64_t;
        using coeff1 = std::int32_t;
        using coeff2 = std::int32_t;

        static constexpr bool k_is_fixed_point   = true;
        static constexpr int  k_state_frac_bits  = 60; ///< Q2.60: two integer guard bits, the top bit spare
        static constexpr int  k_input_shift      = 29; ///< Q0.31 -> Q2.60
        static constexpr int  k_coeff1_frac_bits = 31;
        static constexpr int  k_coeff2_frac_bits = 30;
        static_assert(k_state_frac_bits - k_input_shift == 31, "Q31 ladder: Q0.31 << 29 = Q2.60");

        static constexpr coeff1 make_coeff1(double a) noexcept {
            return detail::round_sat<std::int32_t>(a * static_cast<double>(std::int64_t{1} << k_coeff1_frac_bits));
        }
        static constexpr coeff2 make_coeff2(double c) noexcept {
            return detail::round_sat<std::int32_t>(c * static_cast<double>(std::int64_t{1} << k_coeff2_frac_bits));
        }
        static constexpr state input(std::int32_t x) noexcept {
            return static_cast<std::int64_t>(x) * (std::int64_t{1} << k_input_shift);
        }

        /// (sum of P states) / (P 2^29), round half up, saturating to Q31: the
        /// sum taken as P sum(v >> lg P) + sum(v & (P - 1)), so it never
        /// overflows.
        template <int Lg>
        static constexpr std::int32_t output_mean(std::span<const state> v) noexcept {
            std::int64_t high = 0;
            std::int64_t low  = 0;
            for (const state s : v) {
                high += s >> Lg;
                low += s & ((std::int64_t{1} << Lg) - 1);
            }
            return detail::clamp_sat<std::int32_t>((high + (std::int64_t{1} << 28) + (low >> Lg)) >> 29);
        }
        static constexpr std::int32_t output(state y) noexcept {
            const std::array<state, 1> v{y};
            return output_mean<0>(v);
        }
        static constexpr std::int32_t output_mean2(state a, state b) noexcept {
            const std::array<state, 2> v{a, b};
            return output_mean<1>(v);
        }
        static constexpr std::int32_t output_mean4(state a, state b, state c, state d) noexcept {
            const std::array<state, 4> v{a, b, c, d};
            return output_mean<2>(v);
        }

        /// y = x1 + floor(a (x - y1) / 2^31), the 96-bit product exact.
        static constexpr state allpass1(state x, state x1, state y1, coeff1 a) noexcept {
            return detail::wrap_add(x1, detail::mul_q31_floor(a, x - y1));
        }

        /// y = trunc((x2 2^30 + c2 (x - y2) + c1 (x1 - y1)) / 2^30), exact, toward zero.
        static constexpr state allpass2(state x, state x1, state x2, state y1, state y2, coeff2 c1,
                                        coeff2 c2) noexcept {
            detail::wide_accumulator acc;
            acc.add_shifted30(x2);
            acc.mac(c2, x - y2);
            acc.mac(c1, x1 - y1);
            return acc.shifted<true>(30);
        }

        /// Direct form I, exact, toward zero (as Q15's).
        static constexpr state biquad(state x, state x1, state x2, state y1, state y2,
                                      const std::array<coeff2, 5>& c) noexcept {
            detail::wide_accumulator acc;
            acc.mac(c[0], x);
            acc.mac(c[1], x1);
            acc.mac(c[2], x2);
            acc.mac(c[3], y1);
            acc.mac(c[4], y2);
            return acc.shifted<true>(30);
        }
    };
    // ANCHOR_END: iir_q31

    template <sample_type S>
    using iir_state = typename iir_traits<S>::state;

    // ------------------------------------------------------------------------
    // Path primitives (unchecked: the structure kernels state the precondition)
    // ------------------------------------------------------------------------

    /// y = A(x) through n one-multiplier sections; s holds the n + 1 joint
    /// values one path sample ago.
    template <sample_type S>
    inline iir_state<S> allpass1_path(const typename iir_traits<S>::coeff1* TAP_DSP_RESTRICT a, std::size_t n,
                                      iir_state<S>* TAP_DSP_RESTRICT s, iir_state<S> x) noexcept {
        using tr = iir_traits<S>;
        for (std::size_t i = 0; i < n; ++i) {
            const iir_state<S> y = tr::allpass1(x, s[i], s[i + 1], a[i]);
            s[i]                 = x;
            x                    = y;
        }
        s[n] = x;
        return x;
    }

    /// y = A1(x) through n2 second-order sections (c[2i], c[2i + 1]) = (c1,
    /// c2), then, when n1 == 1, one one-multiplier section a; s1 and s2 hold
    /// the n2 + n1 + 1 joint values one and two path samples ago.
    template <sample_type S>
    inline iir_state<S> allpass2_path(const typename iir_traits<S>::coeff2* TAP_DSP_RESTRICT c, std::size_t n2,
                                      const typename iir_traits<S>::coeff1* TAP_DSP_RESTRICT a, std::size_t n1,
                                      iir_state<S>* TAP_DSP_RESTRICT s1, iir_state<S>* TAP_DSP_RESTRICT s2,
                                      iir_state<S> x) noexcept {
        using tr = iir_traits<S>;
        for (std::size_t i = 0; i < n2; ++i) {
            const iir_state<S> y = tr::allpass2(x, s1[i], s2[i], s1[i + 1], s2[i + 1], c[2 * i], c[2 * i + 1]);
            s2[i]                = s1[i];
            s1[i]                = x;
            x                    = y;
        }
        if (n1 != 0) {
            const iir_state<S> y = tr::allpass1(x, s1[n2], s1[n2 + 1], a[0]);
            s2[n2]               = s1[n2];
            s1[n2]               = x;
            x                    = y;
        }
        s2[n2 + n1] = s1[n2 + n1];
        s1[n2 + n1] = x;
        return x;
    }

    // ------------------------------------------------------------------------
    // Coefficient tables in format S (constexpr, fixed capacity)
    // ------------------------------------------------------------------------

    // ANCHOR: iir_halfband_table
    /// A nonlinear two-path half-band in format S: path p's sections, and the
    /// state a channel needs, (n0 + 1) + (n1 + 1) values (path 0's then
    /// path 1's).
    template <sample_type S>
    struct nonlinear_halfband_table {
        using coeff1 = typename iir_traits<S>::coeff1;

        std::array<coeff1, k_max_halfband_coefficients> path0{};
        std::array<coeff1, k_max_halfband_coefficients> path1{};
        std::size_t                                     size0 = 0;
        std::size_t                                     size1 = 0;

        constexpr nonlinear_halfband_table() = default;
        /// From a design (iir_design.h), each coefficient through make_coeff1.
        constexpr explicit nonlinear_halfband_table(const nonlinear_halfband_design& d) noexcept
            : size0(d.path_size(0))
            , size1(d.path_size(1)) {
            for (std::size_t i = 0; i < size0; ++i) {
                path0[i] = iir_traits<S>::make_coeff1(d.path_coefficient(0, i));
            }
            for (std::size_t i = 0; i < size1; ++i) {
                path1[i] = iir_traits<S>::make_coeff1(d.path_coefficient(1, i));
            }
        }

        /// Multiplies per decimated output (the sections of both paths).
        constexpr std::size_t multiplies() const noexcept { return size0 + size1; }
        /// State values per channel.
        constexpr std::size_t state_size() const noexcept { return size0 + size1 + 2; }
    };
    // ANCHOR_END: iir_halfband_table

    /// A near-linear two-path half-band in format S: path 0 the delay of
    /// `delay` path samples (the caller's history supplies it), path 1 the
    /// allpass of order N as n2 second-order sections and n1 (0 or 1)
    /// one-multiplier sections; a channel's state is 2 (n2 + n1 + 1) values
    /// (s1, then s2).
    template <sample_type S>
    struct near_linear_halfband_table {
        using coeff1 = typename iir_traits<S>::coeff1;
        using coeff2 = typename iir_traits<S>::coeff2;

        static constexpr std::size_t k_max_second_order = 16;

        std::array<coeff2, 2 * k_max_second_order> second{};
        std::array<coeff1, 1>                      first{};
        std::size_t                                n2    = 0;
        std::size_t                                n1    = 0;
        std::size_t                                delay = 0; ///< K, path samples

        constexpr near_linear_halfband_table() = default;
        constexpr explicit near_linear_halfband_table(const near_linear_halfband_design& d) noexcept
            : n2(static_cast<std::size_t>(d.order / 2))
            , n1(static_cast<std::size_t>(d.order % 2))
            , delay(static_cast<std::size_t>(d.delay)) {
            for (std::size_t i = 0; i < 2 * n2; ++i) {
                second[i] = iir_traits<S>::make_coeff2(d.sections[i]);
            }
            if (n1 != 0) {
                first[0] = iir_traits<S>::make_coeff1(d.sections[2 * n2]);
            }
        }

        constexpr std::size_t multiplies() const noexcept { return 2 * n2 + n1; }
        constexpr std::size_t state_size() const noexcept { return 2 * (n2 + n1 + 1); }
    };

    /// A four-path stage in format S: path k's sections; a channel's state is
    /// sum_k (n_k + 1) values, path 0's first.
    template <sample_type S>
    struct fourpath_table {
        using coeff1 = typename iir_traits<S>::coeff1;

        static constexpr std::size_t k_max_per_path = 8;

        std::array<std::array<coeff1, k_max_per_path>, 4> paths{};
        std::array<std::size_t, 4>                        sizes{};

        constexpr fourpath_table() = default;
        constexpr explicit fourpath_table(const fourpath_design& d) noexcept {
            for (std::size_t k = 0; k < 4; ++k) {
                sizes[k] = d.paths[k].size();
                for (std::size_t i = 0; i < sizes[k]; ++i) {
                    paths[k][i] = iir_traits<S>::make_coeff1(d.paths[k][i]);
                }
            }
        }

        constexpr std::size_t multiplies() const noexcept { return sizes[0] + sizes[1] + sizes[2] + sizes[3]; }
        constexpr std::size_t state_size() const noexcept { return multiplies() + 4; }
    };

    /// A biquad / first-order cascade in format S (the guard band filter's
    /// Chebyshev II): per section {b0, b1, b2, -a1, -a2} through make_coeff2
    /// (a first-order section has b2 = a2 = 0). A channel's state is 2 (n +
    /// 1) values in fixed point (direct form I: the joints' last two values,
    /// s1 then s2) and 2 n in floating point (transposed direct form II: z1
    /// then z2 per section).
    template <sample_type S>
    struct biquad_cascade_table {
        using coeff2 = typename iir_traits<S>::coeff2;

        std::array<std::array<coeff2, 5>, k_max_cheby2_sections> sections{};
        std::size_t                                              size = 0;

        constexpr biquad_cascade_table() = default;
        constexpr explicit biquad_cascade_table(std::span<const iir_section> s) noexcept
            : size(s.size() < k_max_cheby2_sections ? s.size() : k_max_cheby2_sections) {
            for (std::size_t i = 0; i < size; ++i) {
                sections[i] = {iir_traits<S>::make_coeff2(s[i].b0), iir_traits<S>::make_coeff2(s[i].b1),
                               iir_traits<S>::make_coeff2(s[i].b2), iir_traits<S>::make_coeff2(-s[i].a1),
                               iir_traits<S>::make_coeff2(-s[i].a2)};
            }
        }
        constexpr explicit biquad_cascade_table(const cheby2_design& d) noexcept
            : biquad_cascade_table(std::span<const iir_section>(d.sections.data(), d.count)) {}

        /// Multiplies per sample: five per section (a first-order section is
        /// run as a biquad with two zero coefficients).
        constexpr std::size_t multiplies() const noexcept { return 5 * size; }
        constexpr std::size_t state_size() const noexcept {
            return iir_traits<S>::k_is_fixed_point ? 2 * (size + 1) : 2 * size;
        }
    };

    // ------------------------------------------------------------------------
    // The structure kernels: one commutator step, one channel
    // ------------------------------------------------------------------------

    // ANCHOR: iir_halfband_kernels
    /// Decimate by 2, one output: y[k] = 1/2 [A0(x[2k]) + A1(x[2k - 1])], the
    /// caller passing x[2k] (even) and x[2k - 1] (odd_prev; zero before the
    /// stream). state: t.state_size() values, zero at reset.
    template <sample_type S>
    inline S halfband_decimate(const nonlinear_halfband_table<S>& t, iir_state<S>* TAP_DSP_RESTRICT state, S even,
                               S odd_prev) noexcept {
        TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED(S);
        using tr              = iir_traits<S>;
        const iir_state<S> p0 = allpass1_path<S>(t.path0.data(), t.size0, state, tr::input(even));
        const iir_state<S> p1 = allpass1_path<S>(t.path1.data(), t.size1, state + t.size0 + 1, tr::input(odd_prev));
        return tr::output_mean2(p0, p1);
    }

    /// Interpolate by 2, one input: out[0] = A0(x) (the even output),
    /// out[1] = A1(x) (the odd one), the interpolator's gain of 2 included.
    template <sample_type S>
    inline void halfband_interpolate(const nonlinear_halfband_table<S>& t, iir_state<S>* TAP_DSP_RESTRICT state, S x,
                                     S* TAP_DSP_RESTRICT out) noexcept {
        TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED(S);
        using tr              = iir_traits<S>;
        const iir_state<S> in = tr::input(x);
        out[0]                = tr::output(allpass1_path<S>(t.path0.data(), t.size0, state, in));
        out[1]                = tr::output(allpass1_path<S>(t.path1.data(), t.size1, state + t.size0 + 1, in));
    }
    // ANCHOR_END: iir_halfband_kernels

    /// The near-linear decimator, one output: y[k] = 1/2 [x[2(k - K)] +
    /// A1(x[2k - 1])], the caller passing the delayed even input from its
    /// history and the odd one.
    template <sample_type S>
    inline S near_linear_decimate(const near_linear_halfband_table<S>& t, iir_state<S>* TAP_DSP_RESTRICT state,
                                  S delayed_even, S odd_prev) noexcept {
        TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED(S);
        using tr             = iir_traits<S>;
        const std::size_t  j = t.n2 + t.n1 + 1;
        const iir_state<S> p1 =
            allpass2_path<S>(t.second.data(), t.n2, t.first.data(), t.n1, state, state + j, tr::input(odd_prev));
        return tr::output_mean2(tr::input(delayed_even), p1);
    }

    /// The near-linear interpolator, one input: out[0] = x[n - K] (the delay
    /// path, a copy of the caller's delayed input), out[1] = A1(x[n]).
    template <sample_type S>
    inline void near_linear_interpolate(const near_linear_halfband_table<S>& t, iir_state<S>* TAP_DSP_RESTRICT state,
                                        S delayed, S x, S* TAP_DSP_RESTRICT out) noexcept {
        TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED(S);
        using tr            = iir_traits<S>;
        const std::size_t j = t.n2 + t.n1 + 1;
        out[0]              = delayed;
        out[1] =
            tr::output(allpass2_path<S>(t.second.data(), t.n2, t.first.data(), t.n1, state, state + j, tr::input(x)));
    }

    // ANCHOR: iir_fourpath_kernels
    /// Decimate by 4, one output: y[n] = 1/4 sum_k A_k(x[4n - k]), the caller
    /// passing x[k] = x[4n - k] for k = 0 .. 3.
    template <sample_type S>
    inline S fourpath_decimate(const fourpath_table<S>& t, iir_state<S>* TAP_DSP_RESTRICT state,
                               const S* TAP_DSP_RESTRICT x) noexcept {
        TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED(S);
        using tr = iir_traits<S>;
        std::array<iir_state<S>, 4> p{};
        std::size_t                 at = 0;
        for (std::size_t k = 0; k < 4; ++k) {
            p[k] = allpass1_path<S>(t.paths[k].data(), t.sizes[k], state + at, tr::input(x[k]));
            at += t.sizes[k] + 1;
        }
        return tr::output_mean4(p[0], p[1], p[2], p[3]);
    }

    /// Interpolate by 4, one input: out[k] = A_k(x) for k = 0 .. 3, the
    /// interpolator's gain of 4 included.
    template <sample_type S>
    inline void fourpath_interpolate(const fourpath_table<S>& t, iir_state<S>* TAP_DSP_RESTRICT state, S x,
                                     S* TAP_DSP_RESTRICT out) noexcept {
        TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED(S);
        using tr              = iir_traits<S>;
        const iir_state<S> in = tr::input(x);
        std::size_t        at = 0;
        for (std::size_t k = 0; k < 4; ++k) {
            out[k] = tr::output(allpass1_path<S>(t.paths[k].data(), t.sizes[k], state + at, in));
            at += t.sizes[k] + 1;
        }
    }
    // ANCHOR_END: iir_fourpath_kernels

    // ANCHOR: iir_biquad_kernel
    /// The biquad / first-order cascade, one sample (1:1): direct form I in
    /// fixed point (one rounding per section, toward zero), transposed direct
    /// form II in floating point.
    template <sample_type S>
    inline S biquad_cascade(const biquad_cascade_table<S>& t, iir_state<S>* TAP_DSP_RESTRICT state, S x) noexcept {
        TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED(S);
        using tr       = iir_traits<S>;
        iir_state<S> v = tr::input(x);
        if constexpr (tr::k_is_fixed_point) {
            iir_state<S>* s1 = state;
            iir_state<S>* s2 = state + t.size + 1;
            for (std::size_t i = 0; i < t.size; ++i) {
                const iir_state<S> y = tr::biquad(v, s1[i], s2[i], s1[i + 1], s2[i + 1], t.sections[i]);
                s2[i]                = s1[i];
                s1[i]                = v;
                v                    = y;
            }
            s2[t.size] = s1[t.size];
            s1[t.size] = v;
        }
        else {
            iir_state<S>* z1 = state;
            iir_state<S>* z2 = state + t.size;
            for (std::size_t i = 0; i < t.size; ++i) {
                const auto&        c = t.sections[i];
                const iir_state<S> y = c[0] * v + z1[i];
                z1[i]                = c[1] * v + c[3] * y + z2[i];
                z2[i]                = c[2] * v + c[4] * y;
                v                    = y;
            }
        }
        return tr::output(v);
    }
    // ANCHOR_END: iir_biquad_kernel

    // ANCHOR: iir_biquad_block_kernel
    namespace detail {

        /// Channels biquad_cascade_block runs together through a block: the
        /// most whose state and coefficients a core's floating registers
        /// hold. Measured (bench/README.md, "The block kernel"): two for
        /// float at one or two sections, one otherwise (float from three
        /// sections, double at any count). The grouping changes no
        /// arithmetic: each channel runs biquad_cascade's expressions on its
        /// own state.
        template <typename S, std::size_t K>
        inline constexpr std::size_t k_biquad_block_group = (std::is_same_v<S, float> && K <= 2) ? 2 : 1;

        /// G channels from channel c0, n interleaved frames, their states and
        /// the coefficients in locals for the block. Out of line on purpose
        /// (biquad_cascade_block's comment).
        template <typename S, std::size_t K, std::size_t Ch, std::size_t G>
        TAP_DSP_IIR_NOINLINE void biquad_cascade_group(const biquad_cascade_table<S>& t, S* state, std::size_t c0,
                                                       const S* in, std::size_t n, S* out) noexcept {
            S cf[K][5];
            for (std::size_t s = 0; s < K; ++s) {
                for (std::size_t j = 0; j < 5; ++j) {
                    cf[s][j] = t.sections[s][j];
                }
            }
            S z[G][2 * K];
            for (std::size_t j = 0; j < G; ++j) {
                for (std::size_t i = 0; i < 2 * K; ++i) {
                    z[j][i] = state[(c0 + j) * 2 * K + i];
                }
            }
            const S* x = in + c0;
            S*       o = out + c0;
            for (std::size_t f = 0; f < n; ++f, x += Ch, o += Ch) {
                for (std::size_t j = 0; j < G; ++j) {
                    S v = x[j];
                    for (std::size_t s = 0; s < K; ++s) {
                        const S y   = cf[s][0] * v + z[j][s];
                        z[j][s]     = cf[s][1] * v + cf[s][3] * y + z[j][K + s];
                        z[j][K + s] = cf[s][2] * v + cf[s][4] * y;
                        v           = y;
                    }
                    o[j] = v;
                }
            }
            for (std::size_t j = 0; j < G; ++j) {
                for (std::size_t i = 0; i < 2 * K; ++i) {
                    state[(c0 + j) * 2 * K + i] = z[j][i];
                }
            }
        }

    } // namespace detail

    /// The floating biquad cascade over a block (1:1): n frames of Ch
    /// interleaved channels through K sections, K and Ch compile-time
    /// constants. Per channel it is biquad_cascade, sample after sample: the
    /// same expressions in the same order on the same state, so with
    /// fp-contraction off (fft.h's policy) the output and the state it leaves
    /// are bit for bit the per-sample kernel's, for any split of the frames
    /// into calls, and the two kernels may be mixed on one state. With
    /// contraction on, a compiler may fuse the two kernels' expressions (and
    /// the single-frame path's, below) differently, and their last bits may
    /// then differ on a target with a fused multiply-add: the contract there
    /// is the tolerance against the double golden model, as for every
    /// floating kernel in this file.
    ///
    ///   - t: the cascade's table; precondition t.size == K (TAP_EXPECTS). In
    ///     a build that does not check, the first K rows of t.sections run
    ///     whatever t.size says: a row past t.size is all zeros, a section
    ///     whose output is zero.
    ///   - state: Ch consecutive per-channel states of 2 K words each, z1[K]
    ///     then z2[K], the layout biquad_cascade reads for t.size == K
    ///     (t.state_size() words per channel); zero-initialized by the
    ///     caller.
    ///   - in, out: n * Ch samples each; the same buffer (in place) or
    ///     disjoint buffers. Neither overlaps state or t.
    ///
    /// What it is for: biquad_cascade called per sample reloads every
    /// coefficient and state word from memory on every sample (a store to
    /// the output may alias both), so a stage that filters a block pays for
    /// memory traffic the recursion does not need. Here a group of channels
    /// (detail::k_biquad_block_group) runs through the whole block with the
    /// coefficients and its states copied into locals, and the states are
    /// written back once. Measured against the per-sample kernel in the
    /// Butterworth stage (bench/README.md, "The block kernel"; float, calls
    /// of 128 frames, 1 to 12 channels, the Cortex-M4F, M33, M55 and
    /// Hexagon): 0.26 to 0.52 of its instructions at orders 2 to 8, 0.58 to
    /// 0.67 at order 16; a few percent on the soft-float Cortex-M4, where
    /// the arithmetic is the cost.
    ///
    /// The group's loop is kept out of line (detail::biquad_cascade_group, a
    /// performance attribute only): inlined into a caller's loop over stages
    /// it took its shape from the caller (on the Cortex-M55 GCC vectorized
    /// across sections there, 1.5 times the instructions of the same stage
    /// called once, and other last bits under contraction). Out of line it
    /// costs a call and a copy of the coefficients per group.
    ///
    /// Honest limit, short calls: the copies are per call and per group, 5 K
    /// coefficients in and 2 K state words in and out per channel, so a call
    /// of a few frames gains less (0.36 to 0.85 at 4 frames), and one frame
    /// through four or more sections would cost more than the per-sample
    /// kernel (measured before this rule: 1.5 times on Hexagon at K = 4, one
    /// channel). A single frame through K >= 4 sections therefore goes to
    /// biquad_cascade itself, at the per-sample kernel's cost within 3 %;
    /// every other call takes the block path, where a single frame costs
    /// 0.63 to 0.84 of the per-sample kernel's at K = 2, but 1.03 on Hexagon
    /// for one channel.
    ///
    /// Floating point only: the fixed-point cascade (direct form I, its state
    /// one word longer per row) has the same reloads and no block kernel yet.
    template <std::floating_point S, std::size_t K, std::size_t Ch>
    inline void biquad_cascade_block(const biquad_cascade_table<S>& t, S* state, const S* in, std::size_t n,
                                     S* out) noexcept {
        static_assert(K >= 1 && K <= k_max_cheby2_sections, "biquad_cascade_block: 1 ... 8 sections");
        static_assert(Ch >= 1, "biquad_cascade_block: at least one channel");
        TAP_DSP_EXPECTS_SUBNORMALS_FLUSHED(S);
        TAP_EXPECTS(t.size == K);
        if constexpr (K >= 4) {
            if (n == 1) {
                for (std::size_t c = 0; c < Ch; ++c) {
                    out[c] = biquad_cascade<S>(t, state + c * 2 * K, in[c]);
                }
                return;
            }
        }
        constexpr std::size_t k_group = detail::k_biquad_block_group<S, K>;
        std::size_t           c0      = 0;
        if constexpr (Ch >= k_group) {
            for (; c0 + k_group <= Ch; c0 += k_group) {
                detail::biquad_cascade_group<S, K, Ch, k_group>(t, state, c0, in, n, out);
            }
        }
        if constexpr (Ch % k_group != 0) {
            detail::biquad_cascade_group<S, K, Ch, Ch % k_group>(t, state, c0, in, n, out);
        }
    }
    // ANCHOR_END: iir_biquad_block_kernel

} // namespace tap::dsp
