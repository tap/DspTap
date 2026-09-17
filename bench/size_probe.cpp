/// @file size_probe.cpp
/// @brief One real-FFT profile alone, for the per-leg .text ceilings.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// docs/audit-fft-and-code-smells.md, Part 10 item 4: a size probe per
// profile, built MinSizeRel like the test legs, measured with
// `arm-none-eabi-size -A` (the .text row alone; Berkeley format folds
// .rodata and NOLOAD sections into "text"). No stdio and no printf, so the
// binary is the startup, the transform and what the transform pulls in;
// the ceilings in .github/workflows/bench.yml are stated against this.
// Same TAP_DSP_SC_PRECISION / TAP_DSP_SC_N selection as icount_main.cpp.
#include <cstddef>

#include "bench_common.h"

#ifndef TAP_DSP_SC_PRECISION
#define TAP_DSP_SC_PRECISION 0
#endif
#ifndef TAP_DSP_SC_N
#define TAP_DSP_SC_N 512
#endif

namespace {

#if TAP_DSP_SC_PRECISION == 0
    using sample = float;
#else
    using sample = double;
#endif

    constexpr std::size_t k_n = TAP_DSP_SC_N;
    sample                g_block[k_n];

} // namespace

int main() {
    tap::dsp::bench::fft_under_test<sample> fft(k_n);
    tap::dsp::bench::xorshift32             rng(0x9E3779B9u);
    rng.fill(g_block, k_n);
    fft.forward_inplace(g_block);
    fft.inverse_inplace(g_block);
    return g_block[0] == g_block[0] ? 0 : 1; // NaN check keeps the transforms live
}
