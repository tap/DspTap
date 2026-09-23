/// @file size_probe.cpp
/// @brief One real-FFT profile alone, for the per-leg .text ceilings.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// docs/audit-fft-and-code-smells.md, Part 10 item 4 (and Part 4's `.text`
// assertion for the float instantiation): a size probe per profile, built
// MinSizeRel like the test legs, measured with `arm-none-eabi-size -A` (the
// .text row alone; Berkeley format folds .rodata and NOLOAD sections into
// "text"). No stdio and no printf, so the binary is the startup, the
// transform and what the transform pulls in; the ceilings in
// .github/workflows/bench.yml are stated against this. Same
// TAP_DSP_SC_PRECISION / TAP_DSP_SC_N selection as icount_main.cpp (0 float,
// 1 double, 2 Q15, 3 Q31).
#include <cstddef>
#include <cstdint>
#include <type_traits>

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
#elif TAP_DSP_SC_PRECISION == 1
    using sample = double;
#elif TAP_DSP_SC_PRECISION == 2
    using sample = std::int16_t;
#elif TAP_DSP_SC_PRECISION == 3
    using sample = std::int32_t;
#else
#error "TAP_DSP_SC_PRECISION must be 0 (float), 1 (double), 2 (Q15) or 3 (Q31)"
#endif

    constexpr std::size_t k_n = TAP_DSP_SC_N;
    sample                g_block[k_n];

    // The transforms' results must reach the exit code, or the linker's
    // --gc-sections could drop what the probe exists to measure. The block is
    // a dependent parameter so each if constexpr branch is instantiated for
    // its own profile only (the float branch's NaN self-compare is a
    // tautology on an integer).
    template <typename Sample>
    int probe(Sample* block) {
        tap::dsp::bench::fft_under_test<Sample> fft(k_n);
        tap::dsp::bench::xorshift32             rng(0x9E3779B9u);
        rng.fill(block, k_n);
        if constexpr (std::is_integral_v<Sample>) {
            const int e = fft.forward_inplace(block) + fft.inverse_inplace(block);
            return (static_cast<int>(block[0]) ^ e) & 1;
        }
        else {
            fft.forward_inplace(block);
            fft.inverse_inplace(block);
            return block[0] == block[0] ? 0 : 1; // NaN check keeps the transforms live
        }
    }

} // namespace

int main() {
    return probe(g_block);
}
