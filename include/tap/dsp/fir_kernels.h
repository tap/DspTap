/// @file fir_kernels.h
/// @brief FIR dot-product kernels: planar, SMLALD dual-MAC, and channel-parallel.
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
// Raspberry Pi Pico 2). Gated off when MVE is present: on M55 the compiler
// already auto-vectorizes the scalar loop with Helium and the intrinsic
// path would replace vectors with dual-MACs (SampleRateTap
// docs/PERFORMANCE.md, hypothesis 4). Bit-exactness: each 16x16 product is
// exact in int32 and the int64 accumulation is associative, so pairing
// changes no output bit.
#if defined(__ARM_FEATURE_DSP) && !defined(__ARM_FEATURE_MVE)
#include <arm_acle.h>
#define TAP_DSP_Q15_SMLALD 1
#else
#define TAP_DSP_Q15_SMLALD 0
#endif
// ANCHOR_END: opt_smlald_gate

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

    // ANCHOR: rs_dot_row
    /// Dot product of a coefficient row against a history window, in the
    /// sample type's accumulator domain (see sample_traits.h): tap-order
    /// accumulation, single rounding in finalize.
    template <sample_type S>
    inline S dot_row(const typename sample_traits<S>::coeff* TAP_DSP_RESTRICT row, const S* TAP_DSP_RESTRICT hist,
                     std::size_t taps) noexcept {
        using tr = sample_traits<S>;
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
            for (; t < taps; ++t) // odd-tap tail
                acc = tr::mac(acc, hist[t], row[t]);
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
