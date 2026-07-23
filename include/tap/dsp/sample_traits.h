/// @file sample_traits.h
/// @brief Sample-format customization point for FIR datapaths: float, Q15, Q31.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Carried from SampleRateTap (include/srt/sample_traits.h), where it
// parameterizes the ASRC datapath; promoted here as the family's sample-format
// substrate, shared with RatioTap's fixed-ratio hot loop and available to any
// primitive growing an embedded fixed-point profile. This header is the
// *format core* only — how each sample type stores coefficients, accumulates
// dot products, and rounds/saturates back to samples. Engine-specific
// extensions (e.g. SampleRateTap's inter-phase coefficient blending) layer on
// top by deriving from these specializations and refining the concept.
//
// Three sample types are provided:
//
//  - float        : float I/O and coefficients, double accumulation
//  - std::int16_t : Q15 samples, Q1.14 coefficients, int64 accumulation,
//                   saturating output
//  - std::int32_t : Q31 samples, Q1.30 coefficients, int64 accumulation,
//                   saturating output
//
// The fixed-point profiles exist for targets where double (and sometimes any
// float) is unaffordable: M33-class cores measured ~19x the instruction count
// running SampleRateTap's float datapath vs Q15 (soft-double accumulation).
// Expected deployments include Bluetooth-adjacent converters (RatioTap) and
// M33/M55-class eurorack/pedal targets (TapTools primitives, opt-in per
// primitive — each adoption is its own documented Q-format design).
//
// Deliberately raw sample types + traits, not wrapper classes with operator
// overloads: the Q-format ladder (headroom bits, accumulator pre-shifts,
// where the single rounding happens) is a per-datapath design decision that
// must stay visible at the use site and pinnable by tests — see the Q31
// mac() below for an example a `q31` class could not express honestly. Raw
// int16_t/int32_t also match what codecs and C ABIs actually deliver, and
// keep memcpy/SIMD access patterns (e.g. the SMLALD pair loads in
// fir_kernels.h) legal.
#pragma once

#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>

namespace tap::dsp {

    namespace detail {

        // ANCHOR: st_roundsat
        /// Round-and-saturate a double to a signed integer coefficient/sample type.
        template <typename I>
        constexpr I round_sat(double v) noexcept {
            constexpr double lo = static_cast<double>(std::numeric_limits<I>::min());
            constexpr double hi = static_cast<double>(std::numeric_limits<I>::max());
            const double     r  = v < 0.0 ? v - 0.5 : v + 0.5; // round half away from zero
            if (r <= lo) {
                return std::numeric_limits<I>::min();
            }
            if (r >= hi) {
                return std::numeric_limits<I>::max();
            }
            return static_cast<I>(r);
        }
        // ANCHOR_END: st_roundsat

        /// Saturate a 64-bit accumulator result to a narrower signed integer.
        template <typename I>
        constexpr I clamp_sat(std::int64_t v) noexcept {
            constexpr auto lo = static_cast<std::int64_t>(std::numeric_limits<I>::min());
            constexpr auto hi = static_cast<std::int64_t>(std::numeric_limits<I>::max());
            return static_cast<I>(v < lo ? lo : (v > hi ? hi : v));
        }

    } // namespace detail

    /// Primary template intentionally undefined; specialize per sample type.
    template <typename T>
    struct sample_traits;

    /// Float datapath: float samples and coefficients, double accumulation.
    /// The double accumulator keeps the dot-product noise floor far below a
    /// 120 dB transparency target; float coefficient storage quantizes the
    /// filter at roughly -150 dB, negligible against the same target.
    template <>
    struct sample_traits<float> {
        using coeff = float;  ///< stored filter coefficient type
        using accum = double; ///< dot-product accumulator type

        /// Convert a double-precision designed coefficient to storage form.
        static coeff make_coeff(double c) noexcept { return static_cast<coeff>(c); }
        /// Coefficient units per 1.0 (used by row-sum-preserving quantization;
        /// unity for floating storage, where no correction runs).
        static constexpr double k_coeff_scale = 1.0;

        /// acc + x * c, in the accumulator domain.
        static accum mac(accum acc, float x, coeff c) noexcept {
            return acc + static_cast<double>(x) * static_cast<double>(c);
        }

        /// Convert the accumulator to an output sample (saturates for fixed point).
        static float finalize(accum acc) noexcept { return static_cast<float>(acc); }

        /// The zero/silence sample value.
        static float silence() noexcept { return 0.0f; }
    };

    // ANCHOR: st_q15_core
    /// Q15 fixed-point datapath (samples are int16_t in Q0.15).
    ///
    /// Coefficients are stored in Q1.14: a unity-DC prototype's peak tap
    /// reaches ~1.0, which does not fit Q0.15, so one headroom bit is traded
    /// for one precision bit. Products are Q0.15 x Q1.14 = Q29 and are summed
    /// exactly in int64 (48-80 taps add ~6-7 bits — no overflow, no
    /// intermediate rounding). The single rounding happens in finalize():
    /// Q29 -> Q15 with round-half-up and saturation. Coefficient quantization
    /// (Q14, ~-86 dB) and output quantization (Q15) set the noise floor — both
    /// at the format's own limit, so a converter built on this is
    /// Q15-transparent.
    template <>
    struct sample_traits<std::int16_t> {
        using coeff = std::int16_t;
        using accum = std::int64_t;

        // ANCHOR: st_q15_coeff
        static coeff make_coeff(double c) noexcept {
            return detail::round_sat<coeff>(c * 16384.0); // Q1.14
        }
        static constexpr double k_coeff_scale = 16384.0; // Q1.14 units per 1.0
        // ANCHOR_END: st_q15_coeff

        // ANCHOR: st_q15_mac
        static accum mac(accum acc, std::int16_t x, coeff c) noexcept {
            return acc + static_cast<std::int64_t>(static_cast<std::int32_t>(x) * static_cast<std::int32_t>(c));
        }
        // ANCHOR_END: st_q15_mac

        // ANCHOR: st_q15_finalize
        static std::int16_t finalize(accum acc) noexcept {
            // Round-half-up, not half-even: the bias is a fraction of one
            // sub-LSB rounding step, far below the Q15 noise floor.
            return detail::clamp_sat<std::int16_t>((acc + (1 << 13)) >> 14); // Q29 -> Q15
        }
        // ANCHOR_END: st_q15_finalize

        static std::int16_t silence() noexcept { return 0; }
    };
    // ANCHOR_END: st_q15_core

    // ANCHOR: st_q31_core
    /// Q31 fixed-point datapath (samples are int32_t in Q0.31).
    ///
    /// Coefficients are stored in Q1.30 (one headroom bit for the ~1.0 peak
    /// tap). A full-precision product would be Q0.31 x Q1.30 = 62 bits, which
    /// overflows int64 once ~48 of them are summed, so each product is
    /// pre-shifted down 16 bits (Q45) before accumulation; the discarded bits
    /// sit 14 bits below the final Q31 LSB, far beneath the format's noise
    /// floor. finalize() rounds Q45 -> Q31 with saturation.
    template <>
    struct sample_traits<std::int32_t> {
        using coeff = std::int32_t;
        using accum = std::int64_t;

        static coeff make_coeff(double c) noexcept {
            return detail::round_sat<coeff>(c * 1073741824.0); // Q1.30
        }
        static constexpr double k_coeff_scale = 1073741824.0; // Q1.30 units per 1.0

        // ANCHOR: st_q31_mac
        static accum mac(accum acc, std::int32_t x, coeff c) noexcept {
            return acc + ((static_cast<std::int64_t>(x) * c) >> 16); // Q61 -> Q45
        }
        // ANCHOR_END: st_q31_mac

        static std::int32_t finalize(accum acc) noexcept {
            return detail::clamp_sat<std::int32_t>((acc + (1 << 13)) >> 14); // Q45 -> Q31
        }

        static std::int32_t silence() noexcept { return 0; }
    };
    // ANCHOR_END: st_q31_core

    // ANCHOR: st_core_concept
    /// Satisfied by any type with a complete format-core sample_traits
    /// specialization — everything the FIR kernels (fir_kernels.h) and table
    /// builders require. Consumers with richer datapaths refine this concept
    /// over their own trait extensions (SampleRateTap's sample_type adds the
    /// inter-phase blend contract on top).
    template <typename T>
    concept sample_type =
        requires(T x, double d, typename sample_traits<T>::accum a, typename sample_traits<T>::coeff c) {
            { sample_traits<T>::make_coeff(d) } -> std::same_as<typename sample_traits<T>::coeff>;
            { sample_traits<T>::mac(a, x, c) } -> std::same_as<typename sample_traits<T>::accum>;
            { sample_traits<T>::finalize(a) } -> std::same_as<T>;
            { sample_traits<T>::silence() } -> std::same_as<T>;
        };

    static_assert(sample_type<float>);
    static_assert(sample_type<std::int16_t>);
    static_assert(sample_type<std::int32_t>);
    // ANCHOR_END: st_core_concept

} // namespace tap::dsp
