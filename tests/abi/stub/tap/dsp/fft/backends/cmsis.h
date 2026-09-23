/// @file cmsis.h
/// @brief TEST STUB for the two-image ABI-tag test: a host-buildable "CMSIS" engine with a different layout.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Shadows include/tap/dsp/fft/backends/cmsis.h for ONE translation unit,
// tests/abi/abi_tag_image.cpp, when tests/CMakeLists.txt builds its second
// image with TAP_DSP_FFT_CMSIS defined and this directory first on the
// include path. It is not CMSIS and computes nothing CMSIS computes: it is
// the split-radix engine behind a deliberately different object layout
// (a 96-byte pad and a call counter), so that the second image's
// basic_real_fft<float> — and everything that embeds one by value — has a
// different layout from the first image's, which is the precondition of
// the cross-image hazard audit F4 describes. Nothing that ships includes
// this file; the real engine's contract numbers are copied so the stub
// satisfies tap::dsp::real_fft_engine exactly as the real one does.

#pragma once

#include <cstddef>

#include "tap/dsp/fft/split_radix.h"

namespace tap::dsp::detail {

    class cmsis_real_fft_f32 {
      public:
        static constexpr std::size_t k_min_size     = 32;
        static constexpr std::size_t k_max_size     = 4096;
        static constexpr bool        k_is_shareable = false;

        explicit cmsis_real_fft_f32(std::size_t n)
            : m_inner(n) {}

        [[nodiscard]] std::size_t size() const noexcept { return m_inner.size(); }

        void forward_inplace(float* a) noexcept {
            m_inner.forward_inplace(a);
            ++m_calls;
        }
        void inverse_inplace(float* a) noexcept {
            m_inner.inverse_inplace(a);
            ++m_calls;
        }

      private:
        [[maybe_unused]] char   m_pad[96]{}; ///< the layout difference, on purpose (unread, by design)
        split_radix_rdft<float> m_inner;
        int                     m_calls = 0;
    };

} // namespace tap::dsp::detail
