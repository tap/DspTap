/// @file fft_arith.h
/// @brief Butterfly arithmetic trait for the real FFT profiles: double, float, Q15 I/O, Q31.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// The sibling of sample_traits.h for the FFT. The FIR trait's mac()/finalize()
// shape (a long exact accumulation with one rounding at the end) is the wrong
// one for a butterfly, which multiplies a sample by a unit-magnitude twiddle
// and adds two samples, stage after stage, with the data staying in its own
// format throughout. This header is the arithmetic the fixed-point kernel
// (fft/fixed_point.h, Stage 3b of docs/audit-fft-and-code-smells.md) is
// written against; the docstrings below are that kernel's specification, and
// tests/test_fft_arith.cpp pins every number in them.
//
// Design (audit Part 7, from Welch 1969 and Oppenheim & Weinstein 1972 on
// finite-register FFT error, and the substrate's rule that a Q ladder is
// visible at the use site):
//
//  - One int32 kernel for both fixed profiles. Q31 data is int32 as is. The
//    Q15 profile widens int16 to int32 on input (<< 14, two guard bits) and
//    narrows on output (round-half-up, saturating, the substrate's finalize
//    rule); fft_arith<std::int16_t> is exactly that pair of I/O-width
//    conversions plus the name of the int32 trait that does the work.
//  - Twiddles are sample_traits<std::int32_t>::coeff, i.e. Q1.30 int32, for
//    BOTH fixed profiles (fft_arith<std::int16_t>::coeff is the same type).
//    1.0 is representable (k_coeff_one = 2^30 < 2^31 - 1), so DC and Nyquist
//    twiddles need no special case; every twiddle is generated in double and
//    rounded once by make_coeff (round_sat: half away from zero, saturating).
//  - The multiply is int32 x Q1.30 -> int64, then >> 30 with ONE round-half-up.
//    That is the single documented rounding point per real product. A
//    complex twiddle multiply is therefore TWO roundings per output
//    component: y_r = sub(mul_coeff(x_r, w_r), mul_coeff(x_i, w_i)),
//    y_i = add(mul_coeff(x_r, w_i), mul_coeff(x_i, w_r)). An
//    accumulate-both-products-in-int64-then-shift-once variant is a
//    different design with different (lower) noise and is NOT bit-identical
//    to this one; the kernel and its Welch-model pins are written against
//    the two-rounding form, and the trait deliberately offers no fused op.
//  - Arm's SMMULR/SMMLAR (32x32, high word, rounded) is NOT a bit-exact seam
//    for mul_coeff: it rounds at bit 32 of the product where mul_coeff rounds
//    at bit 30, and no arrangement reproduces one from the other (shifting
//    the twiddle up two bits makes 1.0 unrepresentable, which is why Q1.30
//    was chosen; shr_round(mul_coeff(x, w), 2) is a double rounding that
//    differs from (x*w + 2^31) >> 32 on about an eighth of all products).
//    An SMMULR build is therefore a separately pinned Arm profile with its
//    own noise numbers, decided only when an on-target measurement justifies
//    it — unlike SMLALD for the FIR kernels, which is exact. No second
//    multiply op is defined here so that the portable arithmetic stays the
//    one contract every host and every QEMU leg runs bit-identically.
//  - Scaling is the kernel's policy; the trait supplies the operations and
//    states the arithmetic conditions under which they suffice:
//      * Fixed scaling is SHIFT-BEFORE-BUTTERFLY: each stage's inputs are
//        shr_round'ed by 2 bits (radix-4), 1 bit (radix-2), and 1 bit before
//        the real post-pass, before the butterfly is computed. Under that
//        ordering the complex magnitude never grows (four inputs each scaled
//        by 1/4, unit-magnitude twiddles), so the magnitude bound is set by
//        the input alone: at most sqrt(2) x full scale for a packed real
//        pair (both components at full scale), plus at most half an LSB per
//        rounding. Shift-AFTER-butterfly is not covered: four rotated
//        full-scale Q2.29 inputs at a 45 degree twiddle sum to 2^31 before
//        the shift and saturate.
//      * Q15 (Q2.29 after widen, full scale 2^29): the bound is 2^29.5, 1.5
//        bits below the int32 rail; the two guard bits suffice with no
//        input-level contract and k_fixed_scaling_input_pre_shift is 0.
//      * Q31 (Q0.31, full scale 2^31 - 1): the bound 2^31.5 exceeds the rail,
//        so fixed scaling takes k_fixed_scaling_input_pre_shift = 1 bit
//        (shr_round, one rounding, -6 dB, one bit of 31) before the first
//        stage, giving a bound of 2^30.5, half a bit below the rail.
//      * Block floating point takes no input pre-shift. Before each stage the
//        kernel reads headroom_bits over the block and right-shifts by
//        max(0, growth bits for the stage - headroom), accumulating the
//        exponent; the trait has no left shift, so BFP never normalizes a
//        quiet block upward, and a shift of 0 with headroom to spare simply
//        keeps the extra precision. The kernel must never let a stage
//        consume the full reported headroom: a block value of -2^20 reports
//        11, and -2^20 << 11 is INT32_MIN, at which sub(0, x) saturates and
//        mul_coeff(x, -1.0) is the one saturating product. Leave at least
//        one bit unused. headroom_bits is 31 for an all-zero block AND for
//        {-1}: 31 means "no information", not "silence".
//  - Twiddles are generated at construction in double, w_k = cos/sin of
//    2*pi*k/N through std::cos/std::sin for the first octant and by exact
//    symmetry for the rest (fft/tables.h), each rounded once by make_coeff.
//    Host libm last-bit differences (glibc, newlib, UCRT, Apple) and
//    fp-contraction of a generator expression (Decision D9; the generators
//    are written so that none is contractible) can move a double that lies
//    within 2^-31 of a Q1.30 rounding boundary onto the other side, so
//    fixed-point outputs are host-identical only if the table is; the 3b
//    battery pins the table's checksum for each certified N, and the
//    transforms' output fingerprints on top, so such a difference is
//    detected rather than silently absorbed. (The quantization bound
//    |w_q - w| <= 0.5 LSB holds on every host.)
//
// Every operation is constexpr and noexcept. The int32 operations preserve
// the data's Q format whatever it is (Q0.31 for the Q31 profile, Q2.29 for
// the widened Q15 profile): only the twiddle's format (Q1.30) enters the
// arithmetic. The floating specializations are the same names over plain
// float/double arithmetic so a profile-generic wrapper can be written once.
#pragma once

#include <bit>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "tap/dsp/sample_traits.h"

namespace tap::dsp {

    /// Primary template intentionally undefined; specialize per sample type.
    template <typename Sample>
    struct fft_arith;

    /// Floating profiles (double is the golden model, float the embedded
    /// profile): every operation is the plain IEEE one, nothing rounds
    /// beyond the format, nothing saturates, and there are no guard bits.
    /// wide == sample, widen/narrow are identities, headroom_bits is 0
    /// because a floating block never needs block scaling.
    template <std::floating_point F>
    struct fft_arith<F> {
        using sample = F;
        using wide   = F; ///< the kernel's data width: the sample itself
        using coeff  = F; ///< twiddle type
        using work   = fft_arith<F>;

        static constexpr bool k_is_fixed_point                = false;
        static constexpr int  k_guard_bits                    = 0;
        static constexpr int  k_coeff_frac_bits               = 0;
        static constexpr int  k_fixed_scaling_input_pre_shift = 0; ///< floating data has an exponent
        static constexpr F    k_coeff_one                     = F{1};

        static constexpr coeff make_coeff(double w) noexcept { return static_cast<coeff>(w); }

        static constexpr wide   widen(sample x) noexcept { return x; }
        static constexpr sample narrow(wide x) noexcept { return x; }

        static constexpr F mul_coeff(F x, coeff w) noexcept { return x * w; }
        static constexpr F add(F a, F b) noexcept { return a + b; }
        static constexpr F sub(F a, F b) noexcept { return a - b; }

        /// x * 2^-bits, exact (a power-of-two scale never rounds in binary
        /// floating point short of underflow).
        /// @pre 0 <= bits <= 31
        static constexpr F shr_round(F x, int bits) noexcept {
            assert(bits >= 0 && bits <= 31);
            F s = F{1};
            for (int i = 0; i < bits; ++i) {
                s *= F{0.5};
            }
            return x * s;
        }

        /// Floating data carries its exponent per element: no block headroom.
        static constexpr int headroom_bits(const F* /*x*/, std::size_t /*n*/) noexcept { return 0; }
    };

    // ANCHOR: fa_q31
    /// The int32 kernel's arithmetic: Q31 data, and the Q15 profile's widened
    /// data. Contracts, per operation (format in -> format out, rounding,
    /// saturation), each pinned by the named test in tests/test_fft_arith.cpp:
    ///
    ///  - make_coeff(double) -> coeff  : Q1.30, round half away from zero,
    ///    saturating (sample_traits<std::int32_t>::make_coeff). 1.0 -> 2^30
    ///    exactly. `TwiddleUnityIsRepresentable`.
    ///  - mul_coeff(x, w) -> sample    : Qm.n x Q1.30 -> Qm.n. Exact int64
    ///    product, then >> 30 with round-half-up (add 2^29, arithmetic shift):
    ///    the single rounding point. Saturates on the way back to int32; with
    ///    |w| <= 1.0 the only input that saturates is INT32_MIN x (-1.0)
    ///    (+2^31 -> INT32_MAX). `MulCoeffIsExactOnRepresentableProducts`,
    ///    `MulCoeffRoundsHalfUp`, `MulCoeffSaturatesAtTheRail`.
    ///  - add(a, b), sub(a, b) -> sample : Qm.n +/- Qm.n -> Qm.n, computed in
    ///    int64 and saturated to [INT32_MIN, INT32_MAX]. Never wraps.
    ///    `AddSubSaturateAtTheRails`.
    ///  - shr_round(x, bits) -> sample : x * 2^-bits, round-half-up in the
    ///    discarded bits (add 2^(bits-1), arithmetic shift, in int64 so the
    ///    rounding add cannot overflow). @pre 0 <= bits <= 31 (asserted);
    ///    bits == 0 is the identity and the result always fits.
    ///    `ShrRoundRoundsHalfUp`.
    ///  - headroom_bits(x, n) -> int   : the number of redundant sign bits
    ///    shared by the whole block — the largest s such that every x[i] << s
    ///    still fits in int32 (Arm's CLS, computed as countl_zero of the OR of
    ///    the one's-complement magnitudes x ^ (x >> 31), minus the sign bit).
    ///    31 for an all-zero or empty block and for {-1} (31 is "no
    ///    information", not a silence signal), 30 for {1}, 0 for INT32_MAX or
    ///    INT32_MIN. A BFP kernel never consumes the full amount (file
    ///    header). `HeadroomBitsOnZeroOneAndFullScale`.
    ///  - k_fixed_scaling_input_pre_shift = 1 : under fixed scaling the Q31
    ///    profile shr_rounds its input by one bit before the first stage
    ///    (file header: the bound 2^31.5 of a full-scale packed pair exceeds
    ///    the rail; with the pre-shift it is 2^30.5). BFP does not use it.
    ///  - widen / narrow                : identities (wide == sample).
    template <>
    struct fft_arith<std::int32_t> {
        using sample = std::int32_t;
        using wide   = std::int32_t;                       ///< the kernel's data width
        using coeff  = sample_traits<std::int32_t>::coeff; ///< Q1.30 twiddle
        using work   = fft_arith<std::int32_t>;
        /// Product width: int32 x Q1.30 is exact in int64.
        using product = std::int64_t;

        static constexpr bool k_is_fixed_point = true;
        static constexpr int  k_guard_bits     = 0; ///< Q0.31 data has no spare width
        /// One bit (-6 dB) before the first stage under fixed scaling: the
        /// price of no guard bits. Not applied under block floating point.
        static constexpr int k_fixed_scaling_input_pre_shift = 1;
        /// Q1.30: the twiddle's fraction bits, shared with the FIR Q31 coefficient.
        static constexpr int   k_coeff_frac_bits = sample_traits<std::int32_t>::k_coeff_frac_bits;
        static constexpr coeff k_coeff_one       = coeff{1} << k_coeff_frac_bits; ///< 1.0 in Q1.30
        static_assert(k_coeff_frac_bits == 30 && k_coeff_one == 1073741824, "twiddles are Q1.30; 1.0 is representable");
        static constexpr int k_sample_bits = 32;

        static constexpr coeff make_coeff(double w) noexcept { return sample_traits<std::int32_t>::make_coeff(w); }

        static constexpr wide   widen(sample x) noexcept { return x; }
        static constexpr sample narrow(wide x) noexcept { return x; }

        static constexpr sample mul_coeff(sample x, coeff w) noexcept {
            const product p = static_cast<product>(x) * static_cast<product>(w);
            return detail::clamp_sat<sample>((p + (product{1} << (k_coeff_frac_bits - 1))) >> k_coeff_frac_bits);
        }

        static constexpr sample add(sample a, sample b) noexcept {
            return detail::clamp_sat<sample>(static_cast<std::int64_t>(a) + static_cast<std::int64_t>(b));
        }

        static constexpr sample sub(sample a, sample b) noexcept {
            return detail::clamp_sat<sample>(static_cast<std::int64_t>(a) - static_cast<std::int64_t>(b));
        }

        static constexpr sample shr_round(sample x, int bits) noexcept {
            assert(bits >= 0 && bits <= 31);
            if (bits == 0) {
                return x;
            }
            return static_cast<sample>((static_cast<std::int64_t>(x) + (std::int64_t{1} << (bits - 1))) >> bits);
        }

        static constexpr int headroom_bits(const sample* x, std::size_t n) noexcept {
            std::uint32_t magnitudes = 0;
            for (std::size_t i = 0; i < n; ++i) {
                // One's-complement magnitude: x for x >= 0, ~x for x < 0, so
                // INT32_MIN maps to 0x7fffffff (zero headroom) and -1 to 0.
                magnitudes |= static_cast<std::uint32_t>(x[i] ^ (x[i] >> (k_sample_bits - 1)));
            }
            return std::countl_zero(magnitudes) - 1;
        }
    };
    // ANCHOR_END: fa_q31

    // ANCHOR: fa_q15
    /// The Q15 profile's I/O width. The kernel does not compute in int16: it
    /// widens each Q0.15 sample to int32 on the way in and narrows on the way
    /// out, and every butterfly runs through fft_arith<std::int32_t> (named
    /// here as `work`). Contracts:
    ///
    ///  - widen(x) -> wide   : Q0.15 -> Q2.29, x << 14. Exact. The two guard
    ///    bits above full scale are what let worst-case radix-4 growth run
    ///    under fixed scaling without an input-level contract — ON THE
    ///    CONDITION that the kernel shifts before each butterfly (2 bits per
    ///    radix-4 stage, 1 per radix-2, 1 before the real post-pass), which
    ///    bounds the complex magnitude by sqrt(2) x 2^29 = 2^29.5, 1.5 bits
    ///    below the int32 rail (file header). No input pre-shift:
    ///    k_fixed_scaling_input_pre_shift = 0. INT16_MIN -> -2^29.
    ///    `WidenPlacesTwoGuardBits`.
    ///  - narrow(y) -> sample : Q2.29 -> Q0.15, (y + 2^13) >> 14 in int64,
    ///    round-half-up, saturating to [INT16_MIN, INT16_MAX] — the
    ///    substrate's finalize rule. `NarrowRoundsHalfUpAndSaturates`.
    ///  - narrow(widen(x)) == x for every int16 value. `WidenNarrowRoundTripsEveryInt16`.
    ///  - coeff, make_coeff, k_coeff_one : the int32 trait's Q1.30 twiddles.
    ///    The Q15 profile does NOT use Q1.14 twiddles.
    template <>
    struct fft_arith<std::int16_t> {
        using sample = std::int16_t;
        using wide   = std::int32_t;                       ///< the kernel's data width
        using coeff  = sample_traits<std::int32_t>::coeff; ///< Q1.30 twiddle, same as Q31
        using work   = fft_arith<std::int32_t>;            ///< the arithmetic the kernel runs on `wide`

        static constexpr bool  k_is_fixed_point                = true;
        static constexpr int   k_guard_bits                    = 2; ///< spare bits above full scale after widen()
        static constexpr int   k_fixed_scaling_input_pre_shift = 0; ///< the guard bits cover the growth bound
        static constexpr int   k_widen_shift                   = 32 - 16 - k_guard_bits; ///< 14: Q0.15 -> Q2.29
        static constexpr int   k_coeff_frac_bits               = work::k_coeff_frac_bits;
        static constexpr coeff k_coeff_one                     = work::k_coeff_one;
        static_assert(k_widen_shift == 14, "Q15 I/O: two guard bits means a 14-bit widen");

        static constexpr coeff make_coeff(double w) noexcept { return work::make_coeff(w); }

        static constexpr wide widen(sample x) noexcept { return static_cast<wide>(x) << k_widen_shift; }

        static constexpr sample narrow(wide x) noexcept {
            return detail::clamp_sat<sample>((static_cast<std::int64_t>(x) + (std::int64_t{1} << (k_widen_shift - 1)))
                                             >> k_widen_shift);
        }
    };
    // ANCHOR_END: fa_q15

    /// Satisfied by the four FFT profiles' sample types: everything a
    /// profile-generic kernel reads — the I/O-width pair, the twiddle type,
    /// its generator and its constants, the scaling constants, and a `work`
    /// trait that carries the butterfly operations over `wide`, all noexcept.
    template <typename T>
    concept fft_sample_type = requires(T x, typename fft_arith<T>::wide y, typename fft_arith<T>::coeff w, double d) {
        typename fft_arith<T>::wide;
        typename fft_arith<T>::coeff;
        typename fft_arith<T>::work;
        requires std::is_same_v<decltype(fft_arith<T>::k_is_fixed_point), const bool>;
        requires std::is_same_v<decltype(fft_arith<T>::k_guard_bits), const int>;
        requires std::is_same_v<decltype(fft_arith<T>::k_coeff_frac_bits), const int>;
        requires std::is_same_v<decltype(fft_arith<T>::k_fixed_scaling_input_pre_shift), const int>;
        requires std::is_same_v<decltype(fft_arith<T>::k_coeff_one), const typename fft_arith<T>::coeff>;
        { fft_arith<T>::make_coeff(d) } noexcept -> std::same_as<typename fft_arith<T>::coeff>;
        { fft_arith<T>::widen(x) } noexcept -> std::same_as<typename fft_arith<T>::wide>;
        { fft_arith<T>::narrow(y) } noexcept -> std::same_as<T>;
        { fft_arith<T>::work::mul_coeff(y, w) } noexcept -> std::same_as<typename fft_arith<T>::wide>;
        { fft_arith<T>::work::add(y, y) } noexcept -> std::same_as<typename fft_arith<T>::wide>;
        { fft_arith<T>::work::sub(y, y) } noexcept -> std::same_as<typename fft_arith<T>::wide>;
        { fft_arith<T>::work::shr_round(y, 1) } noexcept -> std::same_as<typename fft_arith<T>::wide>;
        { fft_arith<T>::work::headroom_bits(&y, std::size_t{1}) } noexcept -> std::same_as<int>;
    };

    static_assert(fft_sample_type<double>);
    static_assert(fft_sample_type<float>);
    static_assert(fft_sample_type<std::int16_t>);
    static_assert(fft_sample_type<std::int32_t>);

    // Both fixed profiles share one twiddle type and 1.0 is representable in it.
    static_assert(std::is_same_v<fft_arith<std::int16_t>::coeff, fft_arith<std::int32_t>::coeff>);
    static_assert(fft_arith<std::int32_t>::make_coeff(1.0) == fft_arith<std::int32_t>::k_coeff_one);
    static_assert(fft_arith<std::int32_t>::k_coeff_one < std::numeric_limits<std::int32_t>::max());

} // namespace tap::dsp
