// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Compile-fail case for chain.h's exact-latency rule (SampleRateTap's
// docs/IIR_PLAN.md, I5): latency_output_frames() compiles only for a chain
// whose every stage reports an exact latency. A recursive stage reports a
// nominal one (nominal_latency_output_frames() and its two siblings), which
// the chain composes into its own nominal, never into the exact fraction.
// tests/CMakeLists.txt builds it through compile_fail/expect_compile_failure.cmake,
// which requires the rule's message as the only error.

#include <cstddef>

#include "tap/dsp/chain.h"
#include "tap/dsp/decimate.h"

namespace {

    class nominal_only_stage {
      public:
        using sample                        = float;
        static constexpr std::size_t k_up   = 1;
        static constexpr std::size_t k_down = 2;
        static constexpr std::size_t channels() noexcept { return 1; }
        std::size_t process(const float* /*in*/, std::size_t n, float* /*out*/) noexcept { return n / 2; }
        std::size_t outputs_for(std::size_t n) const noexcept { return n / 2; }
        void        reset() noexcept {}
        double      nominal_latency_output_frames() const noexcept { return 1.2; }
        double      latency_deviation_output_frames() const noexcept { return 1.28; }
        double      group_delay_output_frames(double /*f*/) const noexcept { return 1.2; }
    };

} // namespace

tap::dsp::exact_ratio exact_latency_of_a_recursive_chain() {
    const tap::dsp::chain<tap::dsp::basic_decimator<float, 2>, nominal_only_stage> c(
        tap::dsp::basic_decimator<float, 2>{}, nominal_only_stage{});
    return c.latency_output_frames();
}
