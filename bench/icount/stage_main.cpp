/// @file stage_main.cpp
/// @brief Instruction-count workloads for the Butterworth stage, gated on their bracketed counts.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// One scenario per binary, selected at compile time (bare metal has no argv):
//
//   TAP_DSP_ST_NAME   the scenario key ("bw4_f32_ch2")
//   TAP_DSP_ST_KIND   0 = streaming: basic_butterworth_lowpass<float, ORDER, CH>
//                         at fc = 0.1 over 16,384 interleaved frames in calls
//                         of TAP_DSP_ST_CALL frames, one bracket per call
//                     1 = set_cutoff(): four redesigns in place, at an
//                         oversampling kernel's cutoffs 0.45 / os for
//                         os = 2, 4, 8, 16, one bracket per call
//                     2 = set_cutoff() and then one read of
//                         latency_deviation_output_frames(), the same four
//                         cutoffs, one bracket per pair: what a caller pays
//                         who redesigns and then asks for the deviation,
//                         wherever the stage computes it
//   TAP_DSP_ST_ORDER  the filter's order (4 or 8)
//   TAP_DSP_ST_CH     channels (kinds 1 and 2: 2)
//   TAP_DSP_ST_CALL   frames per process() call (kind 0)
//
// What the count is: the sum over the brackets (bench/call_marker.h; the
// plugin's TAP_DSP_CALL_MAX total=), not the whole run. The DONE line states
// brackets=<n>, and scripts/icount.py takes the bracket sum as the scenario's
// count and refuses a run whose plugin line is missing or counts another
// number of brackets. A whole-run count would not do: the stage's design at
// construction is soft double on the M4 and the M33, about 17 M instructions
// against 1.8 M of streaming at order 4, so a change of a quarter of the
// kernel's cost would sit inside the ratchet's tolerance. Bracketed, the design, the
// fixture and the checksum fold stay outside, and on the hexagon key so does
// qemu-hexagon's per-process offset (bench/README.md).
//
// Kinds 1 and 2 run the design, which is double: they are gated where the
// target computes double in hardware and built under the informational
// prefix elsewhere (CMakeLists.txt beside this file).
//
// The checksum folds every output word in order (kind 0), or each call's
// window and the bit patterns of its nominal latency and deviation (kinds 1
// and 2), into FNV-1a-64: a change to any output bit or derived value moves
// it. `ok` checks the output is finite and not all zero.
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

#include "bench_common.h"
#include "call_marker.h"
#include "tap/dsp/butterworth.h"

#ifndef TAP_DSP_ST_NAME
#define TAP_DSP_ST_NAME "bw4_f32_ch2"
#endif
#ifndef TAP_DSP_ST_KIND
#define TAP_DSP_ST_KIND 0
#endif
#ifndef TAP_DSP_ST_ORDER
#define TAP_DSP_ST_ORDER 4
#endif
#ifndef TAP_DSP_ST_CH
#define TAP_DSP_ST_CH 2
#endif
#ifndef TAP_DSP_ST_CALL
#define TAP_DSP_ST_CALL 128
#endif

namespace {

    constexpr std::size_t k_ch = TAP_DSP_ST_CH;
    using stage                = tap::dsp::basic_butterworth_lowpass<float, TAP_DSP_ST_ORDER, k_ch>;

    /// One bracketed step: the markers around f().
    template <typename F>
    TAP_DSP_BENCH_NOINLINE void bracketed(F&& f) {
        tap_dsp_call_begin();
        f();
        tap_dsp_call_end();
    }

    template <typename T>
    std::uint64_t fold(std::uint64_t h, T v) {
        unsigned char bytes[sizeof(T)];
        std::memcpy(bytes, &v, sizeof(T));
        for (const unsigned char b : bytes) {
            h ^= b;
            h *= tap::dsp::bench::k_fnv1a64_prime;
        }
        return h;
    }

    struct outcome {
        std::uint64_t checksum;
        bool          ok;
        unsigned      brackets;
    };

#if TAP_DSP_ST_KIND == 0

    constexpr std::size_t k_frames = 16384;
    constexpr std::size_t k_call   = TAP_DSP_ST_CALL;
    static_assert(k_frames % k_call == 0);

    outcome run() {
        stage                       lp(0.1);
        std::vector<float>          in(k_frames * k_ch);
        std::vector<float>          out(k_call * k_ch);
        tap::dsp::bench::xorshift32 rng(0x2545F491u);
        rng.fill(in.data(), in.size());

        std::uint64_t h       = tap::dsp::bench::k_fnv1a64_offset;
        bool          nonzero = false;
        bool          finite  = true;
        for (std::size_t f = 0; f < k_frames; f += k_call) {
            bracketed([&] { lp.process(in.data() + f * k_ch, k_call, out.data()); });
            for (const float v : out) {
                h       = fold(h, v);
                nonzero = nonzero || v != 0.0f;
                finite  = finite && std::isfinite(v);
            }
        }
        return {h, finite && nonzero, static_cast<unsigned>(k_frames / k_call)};
    }

#else

    outcome run() {
        stage         lp(0.45);
        std::uint64_t h      = tap::dsp::bench::k_fnv1a64_offset;
        bool          finite = true;
        unsigned      n      = 0;
        for (int os = 2; os <= 16; os *= 2) {
            const double fc = 0.45 / static_cast<double>(os);
#if TAP_DSP_ST_KIND == 1
            bracketed([&] { lp.set_cutoff(fc); });
            const double dev = lp.latency_deviation_output_frames();
#else
            double dev = 0.0;
            bracketed([&] {
                lp.set_cutoff(fc);
                dev = lp.latency_deviation_output_frames();
            });
#endif
            h      = fold(h, static_cast<std::uint64_t>(lp.window_frames()));
            h      = fold(h, lp.nominal_latency_output_frames());
            h      = fold(h, dev);
            finite = finite && std::isfinite(dev) && lp.window_frames() != 0;
            ++n;
        }
        return {h, finite, n};
    }

#endif

} // namespace

int main() {
    const outcome r = run();
    std::printf("TAP_DSP_ICOUNT_DONE ok=%d engine=butterworth backend=iir scenario=%s brackets=%u checksum=0x%016llx\n",
                r.ok ? 1 : 0, TAP_DSP_ST_NAME, r.brackets, static_cast<unsigned long long>(r.checksum));
    return r.ok ? 0 : 1;
}
