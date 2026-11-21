// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Compile-fail case for the sub-band round trip's whole-delay rule
// (tap/dsp/subband.h; tests/CMakeLists.txt builds it through
// compile_fail/expect_compile_failure.cmake, which requires the build to fail
// with the rule's message as its only error): two stages whose delays are
// compile-time constants (static constexpr latency_output_frames()) and whose
// D tau_down + tau_up is 2 x 1/3 + 0 = 2/3 of a full-rate sample, so the types
// decide that T cannot be exact. Where a stage's delay is decided at run time
// the constructor rejects it instead (test_subband.cpp,
// Subband.NonIntegerDelayIsRejected).

#include <cstddef>

#include "tap/dsp/chain.h"
#include "tap/dsp/subband.h"

namespace {

    // Keeps every second frame; delays a third of an output frame, by type.
    class third_frame_down {
      public:
        using sample                                  = float;
        static constexpr std::size_t           k_up   = 1;
        static constexpr std::size_t           k_down = 2;
        static constexpr std::size_t           channels() noexcept { return 1; }
        static constexpr tap::dsp::exact_ratio latency_output_frames() noexcept { return tap::dsp::exact_ratio{1, 3}; }
        std::size_t                            process(const float* in, std::size_t n, float* out) noexcept {
            std::size_t made = 0;
            for (std::size_t i = 0; i < n; i += 2) {
                out[made++] = in[i];
            }
            return made;
        }
        std::size_t outputs_for(std::size_t n) const noexcept { return (n + 1) / 2; }
        void        reset() noexcept {}
    };

    // Repeats each input twice; no delay.
    class hold_two_up {
      public:
        using sample                                  = float;
        static constexpr std::size_t           k_up   = 2;
        static constexpr std::size_t           k_down = 1;
        static constexpr std::size_t           channels() noexcept { return 1; }
        static constexpr tap::dsp::exact_ratio latency_output_frames() noexcept { return tap::dsp::exact_ratio{0, 1}; }
        std::size_t                            process(const float* in, std::size_t n, float* out) noexcept {
            for (std::size_t i = 0; i < 2 * n; ++i) {
                out[i] = in[i / 2];
            }
            return 2 * n;
        }
        std::size_t outputs_for(std::size_t n) const noexcept { return 2 * n; }
        void        reset() noexcept {}
    };

    class no_change {
      public:
        using sample = float;
        static constexpr std::size_t channels() noexcept { return 1; }
        std::size_t                  latency_frames() const noexcept { return 0; }
        void                         process(const float* /*in*/, std::size_t n, float* change) noexcept {
            for (std::size_t i = 0; i < n; ++i) {
                change[i] = 0.0f;
            }
        }
        void reset() noexcept {}
    };

    // Instantiates the class: the one thing these stages may not do.
    static_assert(tap::dsp::subband_delta<third_frame_down, hold_two_up, no_change>::k_factor == 2);

} // namespace
