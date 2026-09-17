/// @file icount_main.cpp
/// @brief Deterministic fixed workloads for the instruction-count ratchet.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Adapted from MuTap's bench/icount/icount_main.cpp (MIT, MuTap contributors)
// for the real-FFT primitive; the policy is in bench/README.md.
//
// One scenario per binary, selected at compile time because bare-metal
// targets have no argv:
//
//   TAP_DSP_SC_NAME       the scenario key in bench/baselines.json ("rfft_f32_512")
//   TAP_DSP_SC_PRECISION  0 = float (the embedded profile), 1 = double (host-class
//                         targets only: soft-float double is not a profile)
//   TAP_DSP_SC_N          transform size (power of two)
//   TAP_DSP_BENCH_ENGINE  the engine under test (bench_common.h)
//
// The QEMU plugin counts the whole run including construction, so the loop
// is sized so the transforms dominate: k_total_samples samples pass through
// forward + inverse per scenario (2048 iterations at N = 512, 512 at
// N = 2048), against a one-time table build of O(N) trig calls — well under
// 1 % of the total at either size. The input is a small xorshift corpus
// generated once and cycled, so the loop is allocation-free and the only
// non-FFT work per iteration is the copy the class's out-of-place forward()
// and inverse() perform anyway and the checksum fold over every output.
// Nothing in the float scenarios is double: the M4 soft-float leg would
// otherwise measure libgcc.
#include <cinttypes>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

#include "bench_common.h"

#ifndef TAP_DSP_SC_NAME
#define TAP_DSP_SC_NAME "rfft_f32_512"
#endif
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

    constexpr std::size_t k_n             = TAP_DSP_SC_N;
    constexpr std::size_t k_total_samples = std::size_t{1} << 20; // per direction, per scenario
    constexpr std::size_t k_iterations    = k_total_samples / k_n;
    constexpr std::size_t k_corpus_blocks = 4;
    static_assert(k_n >= 4 && (k_n & (k_n - 1)) == 0, "TAP_DSP_SC_N must be a power of two");
    static_assert(k_iterations >= 100, "the loop must dominate construction");

    sample run() {
        tap::dsp::bench::fft_under_test<sample> fft(k_n);

        std::vector<sample>         corpus(k_corpus_blocks * k_n);
        std::vector<sample>         spectrum(k_n);
        std::vector<sample>         out(k_n);
        tap::dsp::bench::xorshift32 rng(0x9E3779B9u);
        rng.fill(corpus.data(), corpus.size());

        sample acc = sample(0);
        for (std::size_t i = 0; i < k_iterations; ++i) {
            const sample* in = corpus.data() + (i % k_corpus_blocks) * k_n;
            fft.forward(in, spectrum.data());
            acc = tap::dsp::bench::fold(acc, spectrum.data(), k_n);
            fft.inverse(spectrum.data(), out.data());
            acc = tap::dsp::bench::fold(acc, out.data(), k_n);
        }
        return acc;
    }

    // The checksum is printed as its bit pattern: exact, and no float-to-double
    // promotion on the way out of a float scenario.
    void print_done(bool ok, sample checksum) {
        if constexpr (sizeof(sample) == sizeof(std::uint32_t)) {
            std::uint32_t bits = 0;
            std::memcpy(&bits, &checksum, sizeof bits);
            std::printf("TAP_DSP_ICOUNT_DONE ok=%d engine=%s scenario=%s checksum=0x%08" PRIx32 "\n", ok ? 1 : 0,
                        tap::dsp::bench::k_engine_name.data(), TAP_DSP_SC_NAME, bits);
        }
        else {
            std::uint64_t bits = 0;
            std::memcpy(&bits, &checksum, sizeof bits);
            std::printf("TAP_DSP_ICOUNT_DONE ok=%d engine=%s scenario=%s checksum=0x%016" PRIx64 "\n", ok ? 1 : 0,
                        tap::dsp::bench::k_engine_name.data(), TAP_DSP_SC_NAME, bits);
        }
    }

} // namespace

int main() {
    const sample checksum = run();
    const bool   ok       = checksum == checksum; // NaN would poison it
    print_done(ok, checksum);
    return ok ? 0 : 1;
}
