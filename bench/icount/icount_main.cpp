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
//   TAP_DSP_SC_PRECISION  0 = float (the embedded floating profile), 1 = double
//                         (host-class targets only: soft-float double is not a
//                         profile), 2 = Q15 (std::int16_t), 3 = Q31
//                         (std::int32_t), the two fixed-point profiles under
//                         scaling::fixed (every leg: fixed point is a profile
//                         on the soft-float M4 too)
//   TAP_DSP_SC_N          transform size (power of two)
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
// the loop is allocation-free. Nothing in the float and fixed-point
// scenarios is double: the M4 soft-float leg would otherwise measure libgcc.
//
// Every output word of every iteration goes through the integer FNV-1a-64
// fold in bench_common.h, so the printed checksum is a bit-exact,
// order-sensitive fingerprint of the engine's output: identical between two
// runs of the same binary, and different for a single 1-ulp change in any
// output. The fixed-point scenarios fold the exponent each transform returns
// as well, so a right bit pattern under a wrong scale would show too. `ok`
// is a sanity check that the last iteration round-trips its input: within
// 1e-3 in the floating scenarios' own precision; in the fixed-point ones
// within eight output LSB referred to the input, i.e. |out * 2^s - x| <=
// 8 * 2^s with s = e_fwd + e_inv + 1 - log2 N, the contract's round-trip
// identity (fft.h; the fixed-point inverse applies no 2/N, the exponents
// carry the scale), evaluated in 64-bit integer arithmetic. Eight is the
// header's pinned worst case for one Q31 fixed transform rounded down to a
// power of two (8.5 LSB at index 0 on the adversarial sweep); this workload
// measures 0.52 (Q15) and 3.6 (Q31) LSB, and a broken transform is off by
// orders of magnitude.
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
#elif TAP_DSP_SC_PRECISION == 1
    using sample = double;
#elif TAP_DSP_SC_PRECISION == 2
    using sample = std::int16_t;
#elif TAP_DSP_SC_PRECISION == 3
    using sample = std::int32_t;
#else
#error "TAP_DSP_SC_PRECISION must be 0 (float), 1 (double), 2 (Q15) or 3 (Q31)"
#endif

    constexpr std::size_t k_n             = TAP_DSP_SC_N;
    constexpr std::size_t k_total_samples = std::size_t{1} << 20; // per direction, per scenario
    constexpr std::size_t k_iterations    = k_total_samples / k_n;
    constexpr std::size_t k_corpus_blocks = 4;
    static_assert(k_n >= 4 && (k_n & (k_n - 1)) == 0, "TAP_DSP_SC_N must be a power of two");
    static_assert(k_iterations >= 100, "the loop must dominate construction");

    struct outcome {
        std::uint64_t checksum;
        bool          round_trips;
    };

// The transform under test is constructed through a function the compiler
// may not inline into run(). The count includes the checksum fold, and the
// fold's register allocation is a function of everything else inlined into
// the same function: the recorded baselines were taken when the CMSIS
// engine's construction reached the workload through an out-of-line
// make_floating_engine<float>, and inlining that construction into run()
// spilled the fold's 64-bit hash and rematerialized the FNV prime on every
// element (+12 % on the m55 float keys with the transform byte-identical;
// bench/README.md, Stage 4). The factory keeps the baseline shape whatever
// the class's constructor does, so the count is the transform plus a fixed
// fold. Measured on the m55 key: -2 / +38 instructions against the recorded
// float baselines and +5 (the call) on each fixed-point scenario.
#if defined(__GNUC__) || defined(__clang__)
#define TAP_DSP_BENCH_NOINLINE [[gnu::noinline]]
#elif defined(_MSC_VER)
#define TAP_DSP_BENCH_NOINLINE __declspec(noinline)
#else
#define TAP_DSP_BENCH_NOINLINE
#endif
    TAP_DSP_BENCH_NOINLINE tap::dsp::bench::fft_under_test<sample> make_fft() {
        return tap::dsp::bench::fft_under_test<sample>(k_n);
    }

// The two workloads are selected by the preprocessor, not by if constexpr
// over one function: the floating run() below is textually the Stage 1b/2b
// workload, and keeping it so is what keeps the recorded float counts at
// +0.00 % (a template with discarded branches was measured to change the
// hot loop's codegen by 4 % on x86-64 and 0.1-0.4 % on the Cortex-M legs).
#if TAP_DSP_SC_PRECISION <= 1

    constexpr sample k_round_trip_tolerance = static_cast<sample>(1e-3);

    outcome run() {
        tap::dsp::bench::fft_under_test<sample> fft = make_fft();

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

#else

    constexpr int log2_of(std::size_t n) noexcept {
        int b = 0;
        while (n > 1) {
            n >>= 1;
            ++b;
        }
        return b;
    }
    constexpr int k_log2_n = log2_of(k_n);

    // The fixed-point workload: both transforms return an exponent, folded
    // after the block it scales; no 2/N anywhere, the exponents carry the
    // scale, and the round-trip check is the contract's identity in int64.
    outcome run() {
        tap::dsp::bench::fft_under_test<sample> fft = make_fft();

        std::vector<sample>         corpus(k_corpus_blocks * k_n);
        std::vector<sample>         spectrum(k_n);
        std::vector<sample>         out(k_n);
        tap::dsp::bench::xorshift32 rng(0x9E3779B9u);
        rng.fill(corpus.data(), corpus.size());

        std::uint64_t h     = tap::dsp::bench::k_fnv1a64_offset;
        const sample* in    = corpus.data();
        int           shift = 0; // e_fwd + e_inv + 1 - log2 N of the last iteration
        for (std::size_t i = 0; i < k_iterations; ++i) {
            in              = corpus.data() + (i % k_corpus_blocks) * k_n;
            const int e_fwd = fft.forward(in, spectrum.data());
            h               = tap::dsp::bench::fold(h, spectrum.data(), k_n);
            h               = tap::dsp::bench::fold(h, &e_fwd, 1);
            const int e_inv = fft.inverse(spectrum.data(), out.data());
            h               = tap::dsp::bench::fold(h, out.data(), k_n);
            h               = tap::dsp::bench::fold(h, &e_inv, 1);
            shift           = e_fwd + e_inv + 1 - k_log2_n;
        }

        bool round_trips = true;
        for (std::size_t i = 0; i < k_n; ++i) {
            // |out * 2^s - x| <= 8 * 2^s: eight output LSB, in the input's units.
            const std::int64_t reconstructed = static_cast<std::int64_t>(out[i]) << shift;
            const std::int64_t err           = reconstructed - static_cast<std::int64_t>(in[i]);
            const std::int64_t tolerance     = std::int64_t{8} << shift;
            if (err > tolerance || err < -tolerance) {
                round_trips = false;
            }
        }
        return {h, round_trips};
    }

#endif

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
