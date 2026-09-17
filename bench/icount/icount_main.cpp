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
// N = 2048), against a one-time table build of O(N) trig calls — under 0.2 %
// of the total at either size on the host (bench/README.md has the
// callgrind figures, including the share of the count that is not the
// transform: the class's out-of-place copies, the 2/N scaling loop and the
// checksum fold, a constant dilution the ratchet's percentages sit on top
// of). The input is a small xorshift corpus generated once and cycled, so
// the loop is allocation-free. Nothing in the float scenarios is double:
// the M4 soft-float leg would otherwise measure libgcc.
//
// Every output word of every iteration goes through the integer FNV-1a-64
// fold in bench_common.h, so the printed checksum is a bit-exact,
// order-sensitive fingerprint of the engine's output: identical between two
// runs of the same binary, and different for a single 1-ulp change in any
// output. `ok` is a sanity check that the last iteration round-trips its
// input, in the scenario's own precision.
#include <cstddef>
#include <cstdint>
#include <cstdio>
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

    constexpr std::size_t k_n                    = TAP_DSP_SC_N;
    constexpr std::size_t k_total_samples        = std::size_t{1} << 20; // per direction, per scenario
    constexpr std::size_t k_iterations           = k_total_samples / k_n;
    constexpr std::size_t k_corpus_blocks        = 4;
    constexpr sample      k_round_trip_tolerance = static_cast<sample>(1e-3);
    static_assert(k_n >= 4 && (k_n & (k_n - 1)) == 0, "TAP_DSP_SC_N must be a power of two");
    static_assert(k_iterations >= 100, "the loop must dominate construction");

    struct outcome {
        std::uint64_t checksum;
        bool          round_trips;
    };

    outcome run() {
        tap::dsp::bench::fft_under_test<sample> fft(k_n);

        std::vector<sample>         corpus(k_corpus_blocks * k_n);
        std::vector<sample>         spectrum(k_n);
        std::vector<sample>         out(k_n);
        tap::dsp::bench::xorshift32 rng(0x9E3779B9u);
        rng.fill(corpus.data(), corpus.size());

        std::uint64_t h  = tap::dsp::bench::k_fnv1a64_offset;
        const sample* in = corpus.data();
        for (std::size_t i = 0; i < k_iterations; ++i) {
            in = corpus.data() + (i % k_corpus_blocks) * k_n;
            fft.forward(in, spectrum.data());
            h = tap::dsp::bench::fold(h, spectrum.data(), k_n);
            fft.inverse(spectrum.data(), out.data());
            h = tap::dsp::bench::fold(h, out.data(), k_n);
        }

        bool round_trips = true;
        for (std::size_t i = 0; i < k_n; ++i) {
            const sample err = out[i] - in[i];
            if (!(err <= k_round_trip_tolerance && err >= -k_round_trip_tolerance)) { // also catches NaN
                round_trips = false;
            }
        }
        return {h, round_trips};
    }

} // namespace

int main() {
    const outcome r = run();
    // %llx with an explicit cast rather than PRIx64: newlib's <inttypes.h>
    // hides the PRI macros from C++ behind __STDC_FORMAT_MACROS.
    std::printf("TAP_DSP_ICOUNT_DONE ok=%d engine=%s backend=%s scenario=%s checksum=0x%016llx\n",
                r.round_trips ? 1 : 0, tap::dsp::bench::k_engine_name, tap::dsp::bench::backend_name<sample>(),
                TAP_DSP_SC_NAME, static_cast<unsigned long long>(r.checksum));
    return r.round_trips ? 0 : 1;
}
