// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Compile-fail case for chain.h's channel rule (tests/CMakeLists.txt builds it
// through compile_fail/expect_compile_failure.cmake, which requires the build
// to fail with the rule's message as its only error): a chain of stages that
// run different channel counts, a mono decimator and a two-channel stage.
// (A stage that declares no count is not a sync_stage at all; test_chain.cpp
// pins that a chain of one is not a type.)

#include <cstddef>

#include "tap/dsp/chain.h"
#include "tap/dsp/decimate.h"

namespace {

    class two_channel_stage {
      public:
        using sample                        = float;
        static constexpr std::size_t k_up   = 1;
        static constexpr std::size_t k_down = 1;
        static constexpr std::size_t channels() noexcept { return 2; }
        std::size_t                  process(const float* /*in*/, std::size_t n, float* /*out*/) noexcept { return n; }
        std::size_t                  outputs_for(std::size_t n) const noexcept { return n; }
        void                         reset() noexcept {}
    };

    static_assert(tap::dsp::chain<tap::dsp::basic_decimator<float, 2>, two_channel_stage>::k_stages == 2);

} // namespace
