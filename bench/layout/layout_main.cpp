/// @file layout_main.cpp
/// @brief Interleaved against planar: the 1:1 chain stages at 1, 2, 6 and 12 channels, counted on every leg.
// SPDX-License-Identifier: MIT
// Copyright 2026 Timothy Place and the DspTap contributors.
//
// Informational workloads, never a gate (bench/README.md, "Layout"). One
// scenario per binary, selected at compile time (bare metal has no argv):
//
//   TAP_DSP_LY_NAME    the scenario key ("chain_f32_ch6_i")
//   TAP_DSP_LY_KIND    0 = basic_butterworth_lowpass<float, 4, Ch> at fc = 0.1
//                      1 = basic_fft_fir<float, Ch>, 127 taps, B = 128 (N = 256)
//                      2 = chain<the Butterworth, the FFT FIR>
//                      3 = a gain stage, one gain per channel (0.5 + 0.01 c
//                          for channel c): stateless, the case where planar
//                          vectorizes over frames and interleaved only over
//                          a channel pattern. Bench-local (bench_gain below),
//                          not a library primitive; planar runs Ch
//                          bench_gain<1>s, each with its own channel's gain.
//   TAP_DSP_LY_CH      channels: 1, 2, 6 (5.1) or 12 (7.1.4)
//   TAP_DSP_LY_PLANAR  0 = one stage of Ch channels over interleaved frames
//                      1 = Ch mono stages, each over its own contiguous
//                          channel buffer (copies of one designed stage, so
//                          both layouts pay one design; the gain stages are
//                          built per channel). Both layouts know their
//                          channel count at compile time.
//
// 16,384 frames per channel in calls of 128 frames, float. Each 128-frame
// step (one process() call interleaved, Ch calls planar) is bracketed by
// bench/call_marker.h's markers, and the figure the job reports is the
// plugin's sum over the brackets (TAP_DSP_CALL_MAX total=), so the designs,
// the fixture and the checksum fold stay outside it, and with them the
// soft-double construction that dominated a whole-run count on the M-class
// keys. The largest bracket (the worst single step: an FFT block lands in
// one) is printed beside it.
//
// The checksum folds each channel's output in frame order into its own
// FNV-1a-64 and then the channel hashes in order, so a pair's two layouts
// print the SAME checksum exactly when their outputs are bit-identical
// (they are wherever the compiler contracts both instantiations alike; GCC's
// default contraction on the M55 does not, for the FFT FIR: bench/README.md,
// "Layout"). `ok` checks the output
// is finite and not all zero. Nothing in the streaming loop is double: the
// M4 soft-float leg would otherwise count libgcc.
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

#include "bench_common.h"
#include "call_marker.h"
#include "tap/dsp/butterworth.h"
#include "tap/dsp/chain.h"
#include "tap/dsp/fft_fir.h"
#include "tap/dsp/kaiser.h"

#ifndef TAP_DSP_LY_NAME
#define TAP_DSP_LY_NAME "chain_f32_ch2_i"
#endif
#ifndef TAP_DSP_LY_KIND
#define TAP_DSP_LY_KIND 2
#endif
#ifndef TAP_DSP_LY_CH
#define TAP_DSP_LY_CH 2
#endif
#ifndef TAP_DSP_LY_PLANAR
#define TAP_DSP_LY_PLANAR 0
#endif

namespace {

    constexpr std::size_t k_ch     = TAP_DSP_LY_CH;
    constexpr std::size_t k_frames = 16384;
    constexpr std::size_t k_call   = 128;
#if TAP_DSP_LY_KIND == 0 || TAP_DSP_LY_KIND == 2
    constexpr double k_fc = 0.1; ///< the Butterworth's cutoff
#endif
#if TAP_DSP_LY_KIND == 1 || TAP_DSP_LY_KIND == 2
    constexpr std::size_t k_taps  = 127; ///< the FFT FIR's taps
    constexpr std::size_t k_block = 128; ///< and its block (N = 256)
#endif

    template <std::size_t Ch>
    using bw = tap::dsp::basic_butterworth_lowpass<float, 4, Ch>;
    template <std::size_t Ch>
    using fir = tap::dsp::basic_fft_fir<float, Ch>;

#if TAP_DSP_LY_KIND == 1 || TAP_DSP_LY_KIND == 2
    std::vector<double> taps() {
        std::vector<double> h(k_taps);
        tap::dsp::design_prototype(h, 1, 0.4, tap::dsp::kaiser_beta(70.0));
        return h;
    }
#endif

    /// A per-channel gain as a 1:1 sync_stage over Ch interleaved channels:
    /// out = in * g[c], g[c] = 0.5 + 0.01 (first + c), the channel count a
    /// template parameter (as every stage's is).
    template <std::size_t Ch>
    class bench_gain {
      public:
        using sample                        = float;
        static constexpr std::size_t k_up   = 1;
        static constexpr std::size_t k_down = 1;
        static constexpr std::size_t channels() noexcept { return Ch; }

        explicit bench_gain(std::size_t first) noexcept {
            for (std::size_t c = 0; c < Ch; ++c) {
                m_g[c] = 0.5f + 0.01f * static_cast<float>(first + c);
            }
        }
        std::size_t process(const float* in, std::size_t n, float* out) noexcept {
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t c = 0; c < Ch; ++c) {
                    out[i * Ch + c] = in[i * Ch + c] * m_g[c];
                }
            }
            return n;
        }
        std::size_t outputs_for(std::size_t n) const noexcept { return n; }
        void        reset() noexcept {}

      private:
        std::array<float, Ch> m_g{};
    };
    static_assert(tap::dsp::sync_stage<bench_gain<12>>);

    template <std::size_t Ch>
    auto make_stage([[maybe_unused]] std::size_t first = 0) {
#if TAP_DSP_LY_KIND == 3
        return bench_gain<Ch>(first);
#elif TAP_DSP_LY_KIND == 0
        return bw<Ch>(k_fc);
#elif TAP_DSP_LY_KIND == 1
        const auto h = taps();
        return fir<Ch>(h, k_block);
#else
        const auto h = taps();
        return tap::dsp::chain<bw<Ch>, fir<Ch>>(bw<Ch>(k_fc), fir<Ch>(h, k_block));
#endif
    }

    /// The input sample of channel c at frame i: xorshift, the same in both layouts.
    float input(std::uint32_t& s) {
        s ^= s << 13;
        s ^= s >> 17;
        s ^= s << 5;
        return static_cast<float>(s) / 4294967296.0f - 0.5f;
    }

    inline std::uint64_t fold1(std::uint64_t h, float v) {
        std::uint32_t bits = 0;
        std::memcpy(&bits, &v, sizeof bits);
        h ^= bits;
        h *= tap::dsp::bench::k_fnv1a64_prime;
        return h;
    }

    /// One bracketed step: the markers around f(). noinline, so callgrind can
    /// collect it on the host (--toggle-collect='*layout_step*') as the
    /// plugin collects the bracket on QEMU.
    template <typename F>
    TAP_DSP_BENCH_NOINLINE void layout_step(F&& f) {
        tap_dsp_call_begin();
        f();
        tap_dsp_call_end();
    }

    struct outcome {
        std::uint64_t checksum;
        bool          ok;
    };

    outcome run() {
        std::array<std::uint64_t, k_ch> h{};
        h.fill(tap::dsp::bench::k_fnv1a64_offset);
        bool          nonzero = false;
        bool          finite  = true;
        std::uint32_t seed    = 0x2545F491u;

#if TAP_DSP_LY_PLANAR == 0
        auto               stage = make_stage<k_ch>();
        std::vector<float> in(k_frames * k_ch), out(k_call * k_ch);
        for (std::size_t i = 0; i < k_frames; ++i) {
            for (std::size_t c = 0; c < k_ch; ++c) {
                in[i * k_ch + c] = input(seed);
            }
        }
        for (std::size_t f = 0; f < k_frames; f += k_call) {
            std::size_t made = 0;
            layout_step([&] { made = stage.process(in.data() + f * k_ch, k_call, out.data()); });
            for (std::size_t i = 0; i < made; ++i) {
                for (std::size_t c = 0; c < k_ch; ++c) {
                    const float v = out[i * k_ch + c];
                    h[c]          = fold1(h[c], v);
                    nonzero       = nonzero || v != 0.0f;
                    finite        = finite && std::isfinite(v);
                }
            }
        }
#else
        const auto                             proto = make_stage<1>();
        std::vector<decltype(make_stage<1>())> stages(k_ch, proto);
#if TAP_DSP_LY_KIND == 3
        for (std::size_t c = 0; c < k_ch; ++c) {
            stages[c] = make_stage<1>(c); // each mono gain stage carries its own channel's gain
        }
#endif
        std::array<std::vector<float>, k_ch> in;
        std::array<std::vector<float>, k_ch> out;
        for (auto& v : in) {
            v.resize(k_frames);
        }
        for (auto& v : out) {
            v.resize(k_call);
        }
        for (std::size_t i = 0; i < k_frames; ++i) {
            for (std::size_t c = 0; c < k_ch; ++c) {
                in[c][i] = input(seed);
            }
        }
        for (std::size_t f = 0; f < k_frames; f += k_call) {
            std::size_t made = 0;
            layout_step([&] {
                for (std::size_t c = 0; c < k_ch; ++c) {
                    made = stages[c].process(in[c].data() + f, k_call, out[c].data());
                }
            });
            for (std::size_t c = 0; c < k_ch; ++c) {
                for (std::size_t i = 0; i < made; ++i) {
                    const float v = out[c][i];
                    h[c]          = fold1(h[c], v);
                    nonzero       = nonzero || v != 0.0f;
                    finite        = finite && std::isfinite(v);
                }
            }
        }
#endif
        std::uint64_t all = tap::dsp::bench::k_fnv1a64_offset;
        for (const std::uint64_t v : h) {
            all ^= v;
            all *= tap::dsp::bench::k_fnv1a64_prime;
        }
        return {all, finite && nonzero};
    }

} // namespace

int main() {
    const outcome r = run();
    std::printf("TAP_DSP_ICOUNT_DONE ok=%d engine=layout backend=%s scenario=%s checksum=0x%016llx\n", r.ok ? 1 : 0,
                tap::dsp::bench::backend_name<float>(), TAP_DSP_LY_NAME, static_cast<unsigned long long>(r.checksum));
    return r.ok ? 0 : 1;
}
