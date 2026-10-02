/// @file fir_kernels.h
/// @brief FIR dot-product kernels: planar, SMLALD dual-MAC, Helium Q15, and channel-parallel.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Carried from SampleRateTap (include/srt/polyphase_filter.h), where these
// kernels are the ASRC's hot loop; promoted here so RatioTap's fixed-ratio
// converter runs the identical, already-measured code. The performance claims
// cited below were measured in SampleRateTap's optimization campaign and are
// regression-gated there by instruction-count CI on Cortex-M33/M55 and
// Hexagon (see that repo's docs/PERFORMANCE.md); a change here lands in every
// consumer on its next submodule bump, so treat the measured comments as
// contracts.
#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>

#include "tap/dsp/sample_traits.h"

// No-alias qualifier for the kernel hot loops: without it the compiler
// versions loops over distinct row/history spans behind a runtime aliasing
// check (verified with -fopt-info-vec; SampleRateTap docs/PERFORMANCE.md,
// hypothesis 2).
#if defined(_MSC_VER)
#define TAP_DSP_RESTRICT __restrict
#else
#define TAP_DSP_RESTRICT __restrict__
#endif

// ANCHOR: opt_smlald_gate
// Dual 16x16 MAC (SMLALD) for the Q15 dot product on Arm cores that have
// the DSP extension but no Helium — the Cortex-M33/M4/M7 class (e.g.
// Raspberry Pi Pico 2). Gated off when MVE is present: the Helium kernel
// below takes the Q15 dot there, eight lanes per reduction where SMLALD
// pairs two (SampleRateTap docs/PERFORMANCE.md, hypothesis 4, measured the
// dual-MAC path against vectors; the next gate records what GCC 13 emits). Bit-exactness: each 16x16 product is
// exact in int32 and the int64 accumulation is associative, so pairing
// changes no output bit.
#if defined(__ARM_FEATURE_DSP) && !defined(__ARM_FEATURE_MVE)
#include <arm_acle.h>
#define TAP_DSP_Q15_SMLALD 1
#else
#define TAP_DSP_Q15_SMLALD 0
#endif
// ANCHOR_END: opt_smlald_gate

// ANCHOR: opt_mve_gate
// Helium (MVE) for the Q15 dot product on Cortex-M55/M85-class cores. The
// claim the SMLALD gate above once relied on — that the compiler
// auto-vectorizes the scalar loop with Helium — does not hold for this loop
// under arm-none-eabi-gcc 13: the int16 x int16 -> int64 reduction compiles
// to one scalar SMLALBB per tap (four instructions per MAC with the loads and
// the low-overhead branch). VMLALDAVA reduces eight 16x16 products into the
// 64-bit accumulator per instruction. Bit-exactness: each product is exact
// and the int64 sum is associative, so lane order changes no output bit (the
// SMLALD argument again). Q15 only: Q31's mac() shifts each product before it
// accumulates, which no Helium reduction reproduces.
#if defined(__ARM_FEATURE_MVE) && (__ARM_FEATURE_MVE & 1)
#include <arm_mve.h>
#define TAP_DSP_Q15_MVE 1
#else
#define TAP_DSP_Q15_MVE 0
#endif
// ANCHOR_END: opt_mve_gate

// Channel-parallel dot product for high channel counts (SampleRateTap
// docs/PERFORMANCE.md, hypothesis C6): history stored frame-major so the
// per-tap inner loop runs across channels — contiguous loads, one
// accumulator lane per channel, coefficient broadcast. Bit-exact because
// each channel's accumulation order over taps is unchanged (lanes are
// channels, not taps), which is what lets the FLOAT path vectorize at all:
// its strict per-channel double accumulation forbids tap-axis SIMD
// (hypothesis 5), but the channel axis is free. Float-only by measurement:
// fixed-point planar dots already auto-vectorize over taps on hosts
// (integer reduction is exactly reassociable) and measured ~1.5x FASTER
// than the channel-parallel form. Host-only: the embedded targets keep
// their proven planar codegen (Helium on M55, SMLALD on M33-class,
// Hexagon's measured scalar floor — hypotheses C4/C5).
#if !defined(__ARM_FEATURE_MVE) && !defined(__ARM_FEATURE_DSP) && !defined(__hexagon__)
#define TAP_DSP_CHANNEL_PARALLEL 1
#else
#define TAP_DSP_CHANNEL_PARALLEL 0
#endif
// Minimum channel count for the frame-major path (overridable for A/B
// measurements; a blend-share planar path stays better at low counts).
#ifndef TAP_DSP_CP_MIN_CHANNELS
#define TAP_DSP_CP_MIN_CHANNELS 4
#endif

namespace tap::dsp {

#if TAP_DSP_Q15_MVE
    namespace detail {
        // ANCHOR: opt_mve_q15
        /// acc + sum_t hist[t] * row[t] over taps Q15 products: eight lanes
        /// per VMLALDAVA over the whole vectors, then one VCTP16-predicated
        /// step (zeroing loads) for the last taps mod 8: no scalar remainder.
        inline std::int64_t mve_q15_accumulate(std::int64_t acc, const std::int16_t* TAP_DSP_RESTRICT row,
                                               const std::int16_t* TAP_DSP_RESTRICT hist, std::size_t taps) noexcept {
            std::size_t t = 0;
            for (; t + 8 <= taps; t += 8) { // whole vectors: no predicate in the loop body
                acc = vmlaldavaq_s16(acc, vld1q_s16(hist + t), vld1q_s16(row + t));
            }
            if (t < taps) { // the last 1 to 7 taps, the dead lanes zeroed
                const mve_pred16_t p = vctp16q(static_cast<std::uint32_t>(taps - t));
                acc                  = vmlaldavaq_s16(acc, vldrhq_z_s16(hist + t, p), vldrhq_z_s16(row + t, p));
            }
            return acc;
        }

        /// acc + sum_t hist[t] * row[taps - 1 - t]: the row read backward
        /// through a gather with descending offsets, eight taps per step; the
        /// last k < 8 taps gather offsets k - 1 ... 0 from the row's start
        /// under a VCTP16 predicate (the inactive lanes' wrapped offsets are
        /// never loaded), so no address outside the row is touched.
        inline std::int64_t mve_q15_accumulate_reversed(std::int64_t acc, const std::int16_t* TAP_DSP_RESTRICT row,
                                                        const std::int16_t* TAP_DSP_RESTRICT hist,
                                                        std::size_t                          taps) noexcept {
            const uint16x8_t down = vddupq_n_u16(7u, 1); // offsets 7, 6, ..., 0
            std::size_t      t    = 0;
            for (; t + 8 <= taps; t += 8) {
                const int16x8_t h = vld1q_s16(hist + t);
                const int16x8_t r = vldrhq_gather_shifted_offset_s16(row + (taps - 8 - t), down);
                acc               = vmlaldavaq_s16(acc, h, r);
            }
            if (t < taps) {
                const auto         k = static_cast<std::uint32_t>(taps - t);
                const mve_pred16_t p = vctp16q(k);
                const int16x8_t    h = vldrhq_z_s16(hist + t, p);
                const int16x8_t    r = vldrhq_gather_shifted_offset_z_s16(row, vddupq_n_u16(k - 1u, 1), p);
                acc                  = vmlaldavaq_s16(acc, h, r);
            }
            return acc;
        }
        // ANCHOR_END: opt_mve_q15
    } // namespace detail
#endif

    // ANCHOR: rs_dot_row
    /// Dot product of a coefficient row against a history window, in the
    /// sample type's accumulator domain (see sample_traits.h): tap-order
    /// accumulation, single rounding in finalize.
    template <sample_type S>
    inline S dot_row(const typename sample_traits<S>::coeff* TAP_DSP_RESTRICT row, const S* TAP_DSP_RESTRICT hist,
                     std::size_t taps) noexcept {
        using tr = sample_traits<S>;
#if TAP_DSP_Q15_MVE
        if constexpr (std::is_same_v<S, std::int16_t>) {
            return tr::finalize(detail::mve_q15_accumulate(0, row, hist, taps));
        }
#endif
#if TAP_DSP_Q15_SMLALD
        if constexpr (std::is_same_v<S, std::int16_t>) {
            std::int64_t acc = 0;
            std::size_t  t   = 0;
            for (; t + 1 < taps; t += 2) {
                // memcpy keeps the 16-bit pair loads alignment-safe; both
                // compile to a single 32-bit load (little-endian packing
                // matches SMLALD's lo/hi lanes).
                std::uint32_t h;
                std::uint32_t r;
                std::memcpy(&h, hist + t, sizeof h);
                std::memcpy(&r, row + t, sizeof r);
                acc = __smlald(static_cast<int16x2_t>(h), static_cast<int16x2_t>(r), acc);
            }
            for (; t < taps; ++t) { // odd-tap tail
                acc = tr::mac(acc, hist[t], row[t]);
            }
            return tr::finalize(acc);
        }
#endif
        typename tr::accum acc{};
        for (std::size_t t = 0; t < taps; ++t) {
            acc = tr::mac(acc, hist[t], row[t]);
        }
        return tr::finalize(acc);
    }
    // ANCHOR_END: rs_dot_row

    // ANCHOR: rs_accumulate_row
    /// dot_row's accumulation without its finalize: acc plus the tap-order
    /// sum of hist[t] * row[t] in the sample type's accumulator domain, so a
    /// caller can sum several rows under ONE rounding point — the branches
    /// of a polyphase decimator whose structural zeros are skipped by
    /// dotting only its nonzero branches (SampleRateTap's rational engine).
    /// Contract: finalize(accumulate_row(accum{}, row, hist, taps)) is
    /// dot_row(row, hist, taps) bit for bit for every sample type: the Q15
    /// path pairs taps with SMLALD behind the same gate (each product exact
    /// in int32, the int64 sum associative), and the floating paths run the
    /// same tap order from the given partial sum, so chaining rows is the
    /// reference mac chain over their concatenation.
    template <sample_type S>
    inline typename sample_traits<S>::accum accumulate_row(typename sample_traits<S>::accum                         acc,
                                                           const typename sample_traits<S>::coeff* TAP_DSP_RESTRICT row,
                                                           const S* TAP_DSP_RESTRICT hist, std::size_t taps) noexcept {
        using tr = sample_traits<S>;
#if TAP_DSP_Q15_MVE
        if constexpr (std::is_same_v<S, std::int16_t>) {
            return detail::mve_q15_accumulate(acc, row, hist, taps);
        }
#endif
#if TAP_DSP_Q15_SMLALD
        if constexpr (std::is_same_v<S, std::int16_t>) {
            std::size_t t = 0;
            for (; t + 1 < taps; t += 2) {
                std::uint32_t h;
                std::uint32_t r;
                std::memcpy(&h, hist + t, sizeof h);
                std::memcpy(&r, row + t, sizeof r);
                acc = __smlald(static_cast<int16x2_t>(h), static_cast<int16x2_t>(r), acc);
            }
            for (; t < taps; ++t) {
                acc = tr::mac(acc, hist[t], row[t]);
            }
            return acc;
        }
#endif
        for (std::size_t t = 0; t < taps; ++t) {
            acc = tr::mac(acc, hist[t], row[t]);
        }
        return acc;
    }
    // ANCHOR_END: rs_accumulate_row

    // ANCHOR: rs_dot_row_reversed
    /// Dot product with the coefficient row read tap-reversed:
    /// y = sum_t hist[t] * row[taps - 1 - t], history walked forward.
    ///
    /// This is the kernel a linear-phase polyphase table halved by symmetry
    /// needs: branch p of a symmetric prototype is branch L-1-p read in
    /// reverse, so the mirrored branches dot the STORED row backward and the
    /// halving costs storage only. Bit-exactness contract: the t-th product
    /// equals dot_row's t-th product against the materialized mirrored row,
    /// and accumulation runs in the same ascending-t order with the same
    /// single finalize — so outputs are bit-identical to keeping the full
    /// table, for every sample type (float's fixed accumulation order
    /// included). On DSP-extension Arm cores the Q15 path pairs taps with the
    /// swapped-lane dual MAC (SMLALDX): each 16x16 product is exact in int32
    /// and the int64 accumulation is associative, so pairing changes no
    /// output bit (the same argument as dot_row's SMLALD gate above).
    template <sample_type S>
    inline S dot_row_reversed(const typename sample_traits<S>::coeff* TAP_DSP_RESTRICT row,
                              const S* TAP_DSP_RESTRICT hist, std::size_t taps) noexcept {
        using tr = sample_traits<S>;
#if TAP_DSP_Q15_MVE
        if constexpr (std::is_same_v<S, std::int16_t>) {
            return tr::finalize(detail::mve_q15_accumulate_reversed(0, row, hist, taps));
        }
#endif
#if TAP_DSP_Q15_SMLALD
        if constexpr (std::is_same_v<S, std::int16_t>) {
            std::int64_t acc = 0;
            std::size_t  t   = 0;
            for (; t + 1 < taps; t += 2) {
                // One 32-bit load per pair on each side; the row load at
                // taps-2-t packs (row[taps-2-t] lo, row[taps-1-t] hi), and
                // SMLALDX's cross pairing multiplies hist.lo * row.hi +
                // hist.hi * row.lo = hist[t]*row[taps-1-t] +
                // hist[t+1]*row[taps-2-t] — exactly the reversed walk.
                std::uint32_t h;
                std::uint32_t r;
                std::memcpy(&h, hist + t, sizeof h);
                std::memcpy(&r, row + (taps - 2 - t), sizeof r);
                acc = __smlaldx(static_cast<int16x2_t>(h), static_cast<int16x2_t>(r), acc);
            }
            for (; t < taps; ++t) { // odd-tap tail
                acc = tr::mac(acc, hist[t], row[taps - 1 - t]);
            }
            return tr::finalize(acc);
        }
#endif
        typename tr::accum acc{};
        for (std::size_t t = 0; t < taps; ++t) {
            acc = tr::mac(acc, hist[t], row[taps - 1 - t]);
        }
        return tr::finalize(acc);
    }
    // ANCHOR_END: rs_dot_row_reversed

    // ANCHOR: opt_dot_tile
    /// One K-channel tile of the channel-parallel dot (hypothesis C6): K
    /// accumulators live in a constexpr-size local array — registers, not
    /// memory — while the tap loop walks the frame-major window with stride
    /// `stride` samples per frame. K is the register-blocking factor; a naive
    /// channels-inner loop with accumulators in memory measures ~2.8x SLOWER
    /// than planar (each mac round-trips its accumulator through the stack).
    template <sample_type S, std::size_t K>
    inline void dot_tile_frame_major(const typename sample_traits<S>::coeff* TAP_DSP_RESTRICT row,
                                     const S* TAP_DSP_RESTRICT x, std::size_t taps, std::size_t stride,
                                     S* TAP_DSP_RESTRICT out) noexcept {
        using tr = sample_traits<S>;
        typename tr::accum acc[K]{};
        for (std::size_t t = 0; t < taps; ++t) {
            const auto                coeff = row[t];
            const S* TAP_DSP_RESTRICT frame = x + t * stride;
            for (std::size_t k = 0; k < K; ++k) {
                acc[k] = tr::mac(acc[k], frame[k], coeff);
            }
        }
        for (std::size_t k = 0; k < K; ++k) {
            out[k] = tr::finalize(acc[k]);
        }
    }
    // ANCHOR_END: opt_dot_tile

    // ANCHOR: rs_dot_rows_frame_major
    // ANCHOR: opt_dot_rows
    /// Channel-parallel dot products over a frame-major history block: all
    /// channels' outputs for one frame in register-blocked tiles of 8/4/2/1.
    /// Per channel the accumulation order over taps equals dot_row's, so the
    /// outputs are bit-exact vs the planar path for every sample type — float
    /// included, since each channel's double accumulator still sums the taps
    /// in the same order (lanes are channels, not taps).
    template <sample_type S>
    inline void dot_rows_frame_major(const typename sample_traits<S>::coeff* TAP_DSP_RESTRICT row,
                                     const S* TAP_DSP_RESTRICT x, std::size_t taps, std::size_t channels,
                                     S* TAP_DSP_RESTRICT out) noexcept {
        std::size_t c = 0;
        for (; c + 8 <= channels; c += 8) {
            dot_tile_frame_major<S, 8>(row, x + c, taps, channels, out + c);
        }
        if (c + 4 <= channels) {
            dot_tile_frame_major<S, 4>(row, x + c, taps, channels, out + c);
            c += 4;
        }
        if (c + 2 <= channels) {
            dot_tile_frame_major<S, 2>(row, x + c, taps, channels, out + c);
            c += 2;
        }
        if (c < channels) {
            dot_tile_frame_major<S, 1>(row, x + c, taps, channels, out + c);
        }
    }
    // ANCHOR_END: rs_dot_rows_frame_major
    // ANCHOR_END: opt_dot_rows

} // namespace tap::dsp
