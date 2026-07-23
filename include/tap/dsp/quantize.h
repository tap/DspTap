/// @file quantize.h
/// @brief Row-sum-preserving coefficient quantization for polyphase tables.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Carried from SampleRateTap, where this correction ran inline in the
// polyphase bank constructor; promoted here because every polyphase table
// with unity per-phase DC gain wants it — RatioTap's fixed-ratio tables
// exactly as much as the ASRC's interpolation bank.
#pragma once

#include <cmath>
#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>

#include "tap/dsp/sample_traits.h"

namespace tap::dsp {

    // ANCHOR: pw_row_sum
    /// Quantizes one polyphase branch's coefficients (double, in units where
    /// 1.0 is unity gain) to the sample type's coefficient format, preserving
    /// the row's DC sum exactly.
    ///
    /// "Note that for DC, this should get you infinite S/N ratio... for every
    /// phase or fractional delay, the FIR coefficients must add to 1."
    ///  -- R. Bristow-Johnson, music-dsp. In double, a well-designed bank
    /// makes every branch's DC sum identical to machine epsilon (zeros at
    /// k*fs ARE branch-DC uniformity, stated in frequency), but independent
    /// per-tap rounding re-breaks it by several LSB. This distributes each
    /// row's total rounding residual to the taps that were rounded furthest
    /// from it (largest-remainder method): every row then sums to
    /// llround(exact_sum * k_coeff_scale), so a table built row-by-row holds
    /// DC gain within one coefficient LSB across all phases.
    ///
    /// For floating-point coefficient types this is a plain conversion (no
    /// correction runs; k_coeff_scale is 1). Design-time code: allocates a
    /// scratch vector for integer formats, so keep it off the audio path.
    ///
    /// \pre dst.size() == src.size()
    template <sample_type S>
    inline void quantize_row_preserving_sum(std::span<const double>                     src,
                                            std::span<typename sample_traits<S>::coeff> dst) {
        using tr            = sample_traits<S>;
        using coeff         = typename tr::coeff;
        const std::size_t n = src.size();
        if constexpr (std::is_floating_point_v<coeff>) {
            for (std::size_t t = 0; t < n; ++t) {
                dst[t] = tr::make_coeff(src[t]);
            }
        }
        else {
            std::vector<double> remainder(n);
            double              exact_sum = 0.0;
            std::int64_t        quant_sum = 0;
            for (std::size_t t = 0; t < n; ++t) {
                const double scaled = src[t] * tr::k_coeff_scale;
                const coeff  q      = tr::make_coeff(src[t]);
                dst[t]              = q;
                remainder[t]        = scaled - static_cast<double>(q);
                exact_sum += scaled;
                quant_sum += static_cast<std::int64_t>(q);
            }
            std::int64_t residual = static_cast<std::int64_t>(std::llround(exact_sum)) - quant_sum;
            while (residual != 0) {
                const double sgn  = residual > 0 ? 1.0 : -1.0;
                std::size_t  best = 0;
                for (std::size_t u = 1; u < n; ++u) {
                    if (sgn * remainder[u] > sgn * remainder[best]) {
                        best = u;
                    }
                }
                dst[best] = static_cast<coeff>(dst[best] + (residual > 0 ? 1 : -1));
                remainder[best] -= sgn;
                residual -= residual > 0 ? 1 : -1;
            }
        }
    }
    // ANCHOR_END: pw_row_sum

} // namespace tap::dsp
